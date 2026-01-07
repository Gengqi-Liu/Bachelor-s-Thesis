function [mEmpfDataBits, data] = EigenModeRx_app(mFrameRxNoCP, params, kanal, len_cInfoBits)
if isfield(params,'V1') && ~isempty(params.V1)
    disp(size(params.V1));
else
    disp('V1 is missing in params.');
end


% EigenModeRx_app
%
% EigenMode receiver (channel estimation + equalization + demapping).
%
% Input:
%   mFrameRxNoCP : [iNfft x iNewNoBlocks x iNoRxAnt] (CP already removed)
%   params       : struct with PHY params; optional V1/U1/S1 for EigenMode
%   kanal        : struct with fields:
%                    schaetzer  ('Zero Forcing' | 'MMSE')
%                    entzerrer  ('Zero Forcing' | 'MMSE' | 'Classical SVD')
%   len_cInfoBits: header bit length (BPSK)
%
% Output:
%   mEmpfDataBits : recovered bitstream (row vector)
%   data          : debug struct (channel estimates, equalized symbols, etc.)

    %% 1) Unpack params
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iModOrd      = params.iModOrd;
    iNoTxAnt     = params.iNoTxAnt;
    iNoRxAnt     = params.iNoRxAnt;
    iNoSubBlocks = params.iNoSubBlocks;

    if isfield(params,'V1') && ~isempty(params.V1)
        V1 = params.V1;  % [Nt x Nt x iNfft]
    else
        warning('EigenModeRx_app: params.V1 missing, using identity precoder (SM fallback).');
        V1 = repmat(eye(iNoTxAnt), 1, 1, iNfft);
    end

    if isfield(params,'U1') && ~isempty(params.U1)
        U1 = params.U1;
    else
        U1 = [];
    end

    if isfield(params,'S1') && ~isempty(params.S1)
        S1 = params.S1;
    else
        S1 = [];
    end

    mFrameRx = mFrameRxNoCP;

    %% 2) Split preamble / data / zero blocks
    AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);

    vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    vPreambleFreq = (1/sqrt(iNfft)) * fft(vPreambleTime, iNfft);

    mPreambleRx = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        mPreambleRx(:,:,:,iCTxAnt) = ...
            mFrameRx(1:iNfft, ...
                     iCTxAnt : (iNoSubBlocks + iNoTxAnt) : ...
                     (iNoSubBlocks*(AnzSubFrames-1) + AnzSubFrames*iNoTxAnt), ...
                     :);
    end

    v = 0;
    mDataRx = zeros(iNfft, iNoBlocks, iNoRxAnt);
    for iSF = 1:AnzSubFrames
        if iSF == AnzSubFrames
            mDataRx(:, iNoSubBlocks*(iSF-1)+1 : iNoSubBlocks*iSF + rem(iNoBlocks,iNoSubBlocks), :) = ...
                mFrameRx(1:iNfft, ...
                         iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v+rem(iNoBlocks,iNoSubBlocks), ...
                         :);
        else
            mDataRx(:, iNoSubBlocks*(iSF-1)+1 : iNoSubBlocks*iSF, :) = ...
                mFrameRx(1:iNfft, ...
                         iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v, ...
                         :);
        end
        v = v + iNoTxAnt + iNoSubBlocks;
    end

    if size(mFrameRx,2) >= 10
        mZeroRx = mFrameRx(:, end-10+1:end, :);
    else
        mZeroRx = mFrameRx;
    end

    %% 3) SNR estimate
    mPreambleRxFreq = (1/sqrt(iNfft)) * fft(mPreambleRx, iNfft, 1);

    mPowNoise    = zeros(iNfft, iNoRxAnt);
    mPowSigNoisy = zeros(iNfft, iNoRxAnt, iNoTxAnt);
    mSNR         = zeros(iNfft, iNoRxAnt, iNoTxAnt);

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            for iSc = 1:iNfft
                mPowNoise(iSc,iCRxAnt) = mean(abs(mZeroRx(iSc,:,iCRxAnt)).^2);
            end
            mPowSigNoisy(:,iCRxAnt,iCTxAnt) = abs(mPreambleRxFreq(:,end,iCRxAnt,iCTxAnt)).^2;
            mSNR(:,iCRxAnt,iCTxAnt) = ...
                (mPowSigNoisy(:,iCRxAnt,iCTxAnt) - mPowNoise(:,iCRxAnt)) ./ ...
                max(mPowNoise(:,iCRxAnt), eps);
        end
    end

    badMask = ~isfinite(mSNR) | (mSNR <= 0);
    mSNR(badMask) = 1e6;

    vSNRperSCperBl = zeros(iNfft,1);
    for iSc = 1:iNfft
        vSNRperSCperBl(iSc) = mean(mean(squeeze(mSNR(iSc,:,:))));
    end

    iSNR = mean(vSNRperSCperBl(:));
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    %% 4) Channel estimation (ZF / MMSE)
    vPreambleFreqCol = vPreambleFreq(:);

    mCTF = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
    mCIR = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);

    if strcmpi(kanal.schaetzer, 'Zero Forcing')
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                for iSF = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iSF,iCRxAnt,iCTxAnt);
                    Hcol  = rxCol ./ vPreambleFreqCol;
                    mCTF(:,iSF,iCRxAnt,iCTxAnt) = Hcol;
                    mCIR(:,iSF,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(Hcol);
                end
            end
        end

    elseif strcmpi(kanal.schaetzer, 'MMSE')
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                snrCol = mSNR(:,iCRxAnt,iCTxAnt);
                den = abs(vPreambleFreqCol).^2 + 1./snrCol;
                den(~isfinite(den) | den<=0) = eps;
                w = conj(vPreambleFreqCol) ./ den;

                for iSF = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iSF,iCRxAnt,iCTxAnt);
                    Hcol  = w .* rxCol;
                    mCTF(:,iSF,iCRxAnt,iCTxAnt) = Hcol;
                    mCIR(:,iSF,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(Hcol);
                end
            end
        end
    else
        error('EigenModeRx_app: unknown channel estimator: %s', kanal.schaetzer);
    end

    %% 5) Per-subcarrier SVD (for debug / Classical SVD fallback)
    CTF = permute(mCTF,[1,4,3,2]);   % [SC x Tx x Rx x SubFrame]
    
    % Use FULL SVD so that U is [Rx x Rx] and fits the preallocation.
    U = zeros(iNoRxAnt, iNoRxAnt, iNfft, AnzSubFrames);
    
    % In full SVD: S is [Rx x Tx] (rectangular). Keep it that way to avoid size mismatch later.
    S = zeros(iNoRxAnt, iNoTxAnt, iNfft, AnzSubFrames);
    
    for j = 1:AnzSubFrames
        for iSc = 1:iNfft
            Htmp = squeeze(CTF(iSc,1:iNoTxAnt,:,j)).';  % [Rx x Tx]
    
            % Full SVD: u is [Rx x Rx], s is [Rx x Tx]
            [u,s,~] = svd(Htmp);
    
            U(:,:,iSc,j) = u;
            S(:,:,iSc,j) = s;
        end
    end

    %% 6) Equalization
    mDataRxFreq = (1/sqrt(iNfft)) * fft(mDataRx, iNfft, 1);  % [SC x Block x Rx]

    if strcmpi(kanal.entzerrer, 'Zero Forcing')
        mDataRxEq = zeros(iNfft, iNoBlocks, iNoTxAnt);

        for iSF = 1:AnzSubFrames
            if iSF == AnzSubFrames
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
            else
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
            end

            for iCB = iCB_range
                for iSc = 1:iNfft
                    H   = squeeze(mCTF(iSc,iSF,:,:));          % [Rx x Tx]
                    Vsc = V1(:,:,iSc);                         % [Tx x Tx]
                    Heff = H * Vsc;                            % [Rx x Tx]
                    y = squeeze(mDataRxFreq(iSc,iCB,:));        % [Rx x 1]
                    mDataRxEq(iSc,iCB,:) = sqrt(iNoTxAnt) * (pinv(Heff) * y);
                end
            end
        end

    elseif strcmpi(kanal.entzerrer, 'MMSE')
        mDataRxEq = zeros(iNfft, iNoBlocks, iNoTxAnt);

        for iSF = 1:AnzSubFrames
            if iSF == AnzSubFrames
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
            else
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
            end

            for iCB = iCB_range
                for iSc = 1:iNfft
                    H   = squeeze(mCTF(iSc,iSF,:,:));
                    Vsc = V1(:,:,iSc);
                    Heff = H * Vsc;
                    y = squeeze(mDataRxFreq(iSc,iCB,:));

                    A = Heff' * Heff + eye(iNoTxAnt) * (1./iSNR);
                    if any(~isfinite(A(:))) || rcond(A) < 1e-10
                        A = Heff' * Heff + eye(iNoTxAnt) * 1e-3;
                    end

                    mDataRxEq(iSc,iCB,:) = sqrt(iNoTxAnt) * (A \ (Heff' * y));
                end
            end
        end

    elseif strcmpi(kanal.entzerrer, 'Classical SVD')
        mDataRxFreqEq = zeros(iNfft, iNoBlocks, iNoRxAnt);
        mDataRxEq     = zeros(iNfft, iNoBlocks, iNoTxAnt);

        if isempty(U1) || isempty(S1)
            U1 = U;
            S1 = S;
            useSubFrameDim = true;
            warning('EigenModeRx_app: params.U1/S1 missing, using per-frame SVD results.');
        else
            useSubFrameDim = (ndims(U1) == 4);
        end

        for b = 1:iNoBlocks
            iSF = ceil(b / iNoSubBlocks);
            if iSF > AnzSubFrames
                iSF = AnzSubFrames;
            end

            for iSc = 1:iNfft
                y = squeeze(mDataRxFreq(iSc,b,:));

                if useSubFrameDim
                    Uuse = U1(:,:,iSc,iSF);
                else
                    Uuse = U1(:,:,iSc);
                end

                mDataRxFreqEq(iSc,b,:) = Uuse' * y;
            end
        end

        for iSF = 1:AnzSubFrames
            if iSF == AnzSubFrames
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
            else
                iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
            end

            for b = iCB_range
                for iSc = 1:iNfft
                    if useSubFrameDim
                        Suse = S1(:,:,iSc,iSF);
                    else
                        Suse = S1(:,:,iSc);
                    end

                    singVals = diag(Suse);
                    singVals(singVals == 0) = eps;

                    z = squeeze(mDataRxFreqEq(iSc,b,1:iNoTxAnt));
                    xEq = sqrt(iNoTxAnt) * (z(:) ./ singVals(1:iNoTxAnt));

                    mDataRxEq(iSc,b,1:iNoTxAnt) = xEq;
                end
            end
        end
    else
        error('EigenModeRx_app: unknown equalizer: %s', kanal.entzerrer);
    end

    %% 7) Symbols -> bits
    vDataRxDet = reshape(mDataRxEq, 1, iNfft*iNoBlocks*iNoTxAnt);
    vDataRxDet(~isfinite(vDataRxDet)) = 0;

    if len_cInfoBits > 0
        vInfoRxDetTp = qamdemod(vDataRxDet(1:len_cInfoBits), 2, 'UnitAveragePower', true);
    else
        vInfoRxDetTp = [];
    end

    if len_cInfoBits < numel(vDataRxDet)
        vDataRxDetTp = qamdemod(vDataRxDet(len_cInfoBits+1:end), 2^iModOrd, 'UnitAveragePower', true);
        mBitsRxDet   = de2bi(vDataRxDetTp, iModOrd, 'right-msb').';
        vBitsRxDet   = mBitsRxDet(:);
    else
        vBitsRxDet = [];
    end

    mEmpfDataBits = [vInfoRxDetTp(:).'  vBitsRxDet(:).'];

    %% 8) Pack debug output
    data = struct();
    data.kanalUeb     = mCTF;
    data.kanalimpuls  = mCIR;
    data.mDataRxEq    = permute(mDataRxEq, [3,1,2]);  % [Tx x SC x Block]
    data.AnzSubFrames = AnzSubFrames;
    data.iSNR         = iSNR;
end