function [mEmpfDataBits, data] = AlamoutiRx_app(mFrameRxNoCP, params, kanal, len_cInfoBits)
%ALAMOUTIRX_APP  Alamouti RX (App version)
% Inputs:
%   mFrameRxNoCP : [iNfft x iNewNoBlocks x iNoRxAnt], CP already removed
%   params       : struct, must contain at least
%                  .iNoBlocks, .iNfft, .iNg, .iNb, .iNoTxAnt, .iNoRxAnt,
%                  .fBBFreq, .fDACFreq, .fCarrFreq
%   kanal        : struct, fields:
%                  .schaetzer : 'Zero Forcing' | 'MMSE'
%   len_cInfoBits: number of BPSK control bits
%
% Outputs:
%   mEmpfDataBits : row vector of recovered bits (header + data)
%   data          : debug struct (channels, equalized symbols, etc.)

    %% 1. Unpack parameters
    iNoBlocks = params.iNoBlocks;
    iNfft     = params.iNfft;
    iNoTxAnt  = params.iNoTxAnt;
    iNoRxAnt  = params.iNoRxAnt;

    fBBFreq   = params.fBBFreq;
    fDACFreq  = params.fDACFreq;
    fCarrFreq = params.fCarrFreq;  %#ok<NASGU>

    if iNoTxAnt ~= 2
        error('AlamoutiRx_app: only 2x? Alamouti supported (iNoTxAnt must be 2).');
    end

    [~, iNewNoBlocks, ~] = size(mFrameRxNoCP);

    mEmpfDataBits = [];
    data          = struct();

    %% 2. Parse Alamouti frame structure (preamble / data / zeros)
    % TX frame (per logical block, nColsPerBlock columns):
    %   generateTx_Alamouti_app uses:
    %     nColsPerBlock = iNoTxAnt + 3 (=5 for 2x?)
    %     vInd = (nColsPerBlock*(iCB-1)+1) : (nColsPerBlock*iCB)
    %     preambles:
    %       Tx1: vInd(1) + 0
    %       Tx2: vInd(2) + 2
    %     data (two Alamouti time slots):
    %       vInd(end)-3 : vInd(end)-2
    %     trailing zero:
    %       vInd(end)
    nColsPerBlock = iNoTxAnt + 3;   % = 5
    neededCols    = nColsPerBlock * iNoBlocks;

    if iNewNoBlocks < neededCols
        error('AlamoutiRx_app: not enough OFDM blocks, need at least %d, got %d.', ...
              neededCols, iNewNoBlocks);
    end

    % preamble: [Nfft x NoBlocks x Nr x Nt]
    mPreambleRx = zeros(iNfft, iNoBlocks, iNoRxAnt, iNoTxAnt);
    % data: [Nfft x (2*NoBlocks) x Nr]
    mDataRx     = zeros(iNfft, 2*iNoBlocks, iNoRxAnt);
    % zeros (noise estimation): [Nfft x NoBlocks x Nr]
    mZeroRx     = zeros(iNfft, iNoBlocks, iNoRxAnt);

    for iCB = 1:iNoBlocks
        vInd = (nColsPerBlock*(iCB-1)+1) : (nColsPerBlock*iCB);

        % preambles
        for iCAnt = 1:iNoTxAnt
            colP = vInd(iCAnt) + (iCAnt-1)*2;
            mPreambleRx(:, iCB, :, iCAnt) = mFrameRxNoCP(:, colP, :);
        end

        % two Alamouti time slots
        mDataRx(:, 2*iCB-1 : 2*iCB, :) = ...
            mFrameRxNoCP(:, vInd(end)-3 : vInd(end)-2, :);

        % trailing zero column
        mZeroRx(:, iCB, :) = mFrameRxNoCP(:, vInd(end), :);
    end

    %% 3. Ideal preamble and SNR estimation
    % Use the exact TX preamble frequency reference if provided
    if isfield(params,'vPreambleFreq') && ~isempty(params.vPreambleFreq)
        vPreambleFreq = params.vPreambleFreq(:);
    else
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
        vPreambleFreq = (1/sqrt(iNfft)) * fft(vPreambleTime(:), iNfft);
    end

    % RX preamble FFT: [Nfft x NoBlocks x Nr x Nt]
    mPreambleRxFreq = 1/sqrt(iNfft) * fft(mPreambleRx, iNfft, 1);

    % noise power: mPowNoise(k, r)
    mPowNoise = zeros(iNfft, iNoRxAnt);
    for iCRxAnt = 1:iNoRxAnt
        z = squeeze(mZeroRx(:,:,iCRxAnt));            % [Nfft x NoBlocks]
        mPowNoise(:,iCRxAnt) = mean(abs(z).^2, 2);    % over blocks
    end

    % signal+noise power: mPowSigNoisy(k, r, t)
    mPowSigNoisy = zeros(iNfft, iNoRxAnt, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            y = squeeze(mPreambleRxFreq(:,:,iCRxAnt,iCTxAnt));    % [Nfft x NoBlocks]
            mPowSigNoisy(:,iCRxAnt,iCTxAnt) = mean(abs(y).^2, 2);
        end
    end

    % SNR per (k,r,t)
    mSNR = zeros(iNfft, iNoRxAnt, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            num = mPowSigNoisy(:,iCRxAnt,iCTxAnt) - mPowNoise(:,iCRxAnt);
            den = mPowNoise(:,iCRxAnt);
            den(den <= 0) = eps;
            snrCol = num ./ den;
            snrCol(~isfinite(snrCol) | snrCol <= 0) = 1e6;
            mSNR(:,iCRxAnt,iCTxAnt) = snrCol;
        end
    end

    iSNR = mean(mSNR(:));
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    %% 4. Channel estimation: ZF / MMSE
    % Outputs: mCTF [Nfft x NoBlocks x Nr x Nt], mCIR same dims
    mCTF = zeros(iNfft, iNoBlocks, iNoRxAnt, iNoTxAnt);
    mCIR = zeros(iNfft, iNoBlocks, iNoRxAnt, iNoTxAnt);

    vPreambleFreqCol = vPreambleFreq(:);  % [Nfft x 1]

    useMMSE = strcmpi(kanal.schaetzer, 'MMSE');

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            snrCol = mSNR(:,iCRxAnt,iCTxAnt);  % [Nfft x 1]

            for iCB = 1:iNoBlocks
                rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);   % [Nfft x 1]

                if ~useMMSE
                    Hcol = rxCol ./ vPreambleFreqCol;
                else
                    den = abs(vPreambleFreqCol).^2 + 1./snrCol;
                    den(den <= 0 | ~isfinite(den)) = eps;
                    w   = conj(vPreambleFreqCol) ./ den;
                    Hcol = w .* rxCol;
                end

                mCTF(:,iCB,iCRxAnt,iCTxAnt) = Hcol;
                mCIR(:,iCB,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(Hcol);
            end
        end
    end

    %% 5. Data FFT
    % mDataRx: [Nfft x 2*NoBlocks x Nr]
    mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);

    %% 6. Build effective Alamouti channel Heff
    % mCTFeff(k,b, 1:2*Nr, 1:2):
    %   upper:     [ H(:,1),  H(:,2)]
    %   lower:     [conj(H(:,2)), -conj(H(:,1))]
    mCTFeff = zeros(iNfft, iNoBlocks, 2*iNoRxAnt, iNoTxAnt);

    for iCB = 1:iNoBlocks
        for k = 1:iNfft
            H = squeeze(mCTF(k,iCB,:,:));  % [Nr x 2]

            h1 = H(:,1);
            h2 = H(:,2);

            % upper part
            mCTFeff(k,iCB,1:iNoRxAnt,1) = h1;
            mCTFeff(k,iCB,1:iNoRxAnt,2) = h2;

            % lower part
            mCTFeff(k,iCB,iNoRxAnt+1:2*iNoRxAnt,1) = conj(h2);
            mCTFeff(k,iCB,iNoRxAnt+1:2*iNoRxAnt,2) = -conj(h1);
        end
    end

    %% 7. Build stacked RX vectors y = [r1; r2*]
    % Conjugate even time-slot symbols, then stack per block
    mDataConjRx              = mDataRxFreq;
    mDataConjRx(:,2:2:end,:) = conj(mDataRxFreq(:,2:2:end,:));

    % mDataNewRx: [Nfft x NoBlocks x 2*Nr]
    mDataNewRx = zeros(iNfft, iNoBlocks, 2*iNoRxAnt);

    for iCB = 1:iNoBlocks
        mDataNewRx(:,iCB,1:iNoRxAnt)           = mDataConjRx(:,2*iCB-1,:);
        mDataNewRx(:,iCB,iNoRxAnt+1:2*iNoRxAnt) = mDataConjRx(:,2*iCB,:);
    end

    %% 8. Alamouti linear detection
    mDataRxEq = zeros(iNoTxAnt, iNfft, iNoBlocks);
    
    for iCB = 1:iNoBlocks
        for k = 1:iNfft
            Heff   = squeeze(mCTFeff(k,iCB,:,:));      % [2*Nr x 2]
            yStack = squeeze(mDataNewRx(k,iCB,:));     % [2*Nr x 1]
    
            % Extract h1 and h2 from the "upper" part of Heff
            h1 = Heff(1:iNoRxAnt, 1);
            h2 = Heff(1:iNoRxAnt, 2);
    
            alpha = sum(abs(h1).^2 + abs(h2).^2);
            if ~isfinite(alpha) || alpha < 1e-12
                alpha = 1e-12;
            end
    
            xHat = (sqrt(iNoTxAnt) / alpha) * (Heff' * yStack);   % [2 x 1]
            mDataRxEq(:,k,iCB) = xHat;
        end
    end


    %% 9. Symbols -> bits (header + payload)
    % Serialize symbols in the same order as SM/VBLAST:
    % [Nfft x NoBlocks x Tx] -> column-major vector
    mDataRxDet = permute(mDataRxEq, [2 3 1]);   % [Nfft x NoBlocks x Tx]
    vDataRxDet = reshape(mDataRxDet, 1, []);    % row vector
    vDataRxDet(~isfinite(vDataRxDet)) = 0;
    vDataRxDet(~isfinite(vDataRxDet)) = 0;

    % header: BPSK
    if len_cInfoBits > 0
        vInfoRxDetTp = qamdemod(vDataRxDet(1:len_cInfoBits), 2, ...
                                'UnitAveragePower', true);
    else
        vInfoRxDetTp = [];
    end

    if isfield(params,'iModOrd')
        iModOrd = params.iModOrd;
    else
        iModOrd = 2;
    end

    % payload: M-QAM
    if len_cInfoBits < numel(vDataRxDet)
        vDataPart    = vDataRxDet(len_cInfoBits+1:end);
        vDataRxDetTp = qamdemod(vDataPart, 2^iModOrd, ...
                                'UnitAveragePower', true);
        mBitsRxDet   = de2bi(vDataRxDetTp, iModOrd, 'right-msb').';
        vBitsRxDet   = mBitsRxDet(:);
    else
        vBitsRxDet = [];
    end

    vBitsAll      = [vInfoRxDetTp(:); vBitsRxDet(:)];
    mEmpfDataBits = vBitsAll.'; 

    %% 10. Fill debug struct
    data.kanalUeb      = mCTF;
    data.kanalimpuls   = mCIR;
    data.mDataRxEq     = mDataRxEq;    % [Tx x Nfft x NoBlocks]

    data.ta_TP         = 1 / fBBFreq;
    data.iNfft         = iNfft;
    data.iModOrd       = iModOrd;
    data.fDACFreq      = fDACFreq;
    data.fBBFreq       = fBBFreq;
    data.fCarrFreq     = fCarrFreq;
    data.iNoBlocks     = iNoBlocks;
    data.iNoTxAnt      = iNoTxAnt;
    data.iNoRxAnt      = iNoRxAnt;
    data.mFrameRxNoCP  = mFrameRxNoCP;
end