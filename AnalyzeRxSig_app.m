function [mEmpfDataBits, data] = AnalyzeRxSig_app(rxFrame, params, kanal)
%ANALYZERXSIG_APP
% RX top-level:
%   1) Downconvert to baseband
%   2) Schmidl & Cox timing
%   3) Cut OFDM blocks and remove CP
%   4) Dispatch to MIMO-mode-specific RX:
%        - AlamoutiRx_app
%        - SM_VBLAST_Rx_app
%        - EigenModeRx_app
%
% Inputs:
%   rxFrame : [Nr x Nsamp] time-domain RX signal
%   params  : physical layer parameter struct
%   kanal   : struct (len_cInfoBits / schaetzer / entzerrer / code / ...)
%
% Outputs:
%   mEmpfDataBits : recovered bits (row vector)
%   data          : debug struct (channels, equalized symbols, ...)

    %% 1. Unpack parameters
    iNoBlocks     = params.iNoBlocks;
    iNfft         = params.iNfft;
    iNg           = params.iNg;
    iNb           = params.iNb;
    iModOrd       = params.iModOrd;
    iNoTxAnt      = params.iNoTxAnt;
    iNoRxAnt      = params.iNoRxAnt;
    iNewNoBlocks  = params.iNewNoBlocks;
    iNoSubBlocks  = params.iNoSubBlocks;

    fBBFreq       = params.fBBFreq;
    fDACFreq      = params.fDACFreq;
    fCarrFreq     = params.fCarrFreq;

    mimoModeStr   = string(params.mimoMode);

    len_cInfoBits = kanal.len_cInfoBits;

    mEmpfDataBits = [];
    data          = struct();

    %% 2. Downconvert to baseband
    % rxFrame: [Nr x Nsamp] at fDACFreq
    mRecFrame = rxFrame;

    vBB1  = DeModulateSignal_app(mRecFrame(1,:), fBBFreq, fDACFreq, fCarrFreq);
    Lbb   = numel(vBB1);

    mDataRxDem = zeros(iNoRxAnt, Lbb);
    mDataRxDem(1,:) = vBB1;

    for iC = 2:iNoRxAnt
        vBB = DeModulateSignal_app(mRecFrame(iC,:), fBBFreq, fDACFreq, fCarrFreq);
        if numel(vBB) ~= Lbb
            L = min(Lbb, numel(vBB));
            mDataRxDem(iC,1:L) = vBB(1:L);
        else
            mDataRxDem(iC,:) = vBB;
        end
    end

    %% 3. Schmidl & Cox timing
    vFrStart   = [];
    vCarOff    = [];
    MetricData = cell(1, iNoRxAnt);

    for iC = 1:iNoRxAnt
        [frameStartTmp, fCarOffTmp, metricStruct] = ...
            EstFrameStart_app(mDataRxDem(iC,:), iNfft, iNg);

        vFrStart = [vFrStart, frameStartTmp];
        vCarOff  = [vCarOff,  fCarOffTmp]; %#ok<AGROW>
        MetricData{iC} = metricStruct;
    end

    if isempty(vFrStart)
        iFrSync = 1;
    else
        iFrSync = min(vFrStart);  % earliest sync point among RX channels
    end

    % S&C training symbol length = iNfft, no CP
    % data block CP start = iFrSync + iNfft
    iFrStart = iFrSync + iNfft;
    iFrStart = max(1, min(iFrStart, Lbb));

    %% 4. Cut OFDM blocks and remove CP
    maxSamplesAvail = Lbb - iFrStart + 1;
    maxBlocksAvail  = floor(maxSamplesAvail / iNb);

    if maxBlocksAvail <= 0
        error(['AnalyzeRxSig_app: not enough samples for one OFDM block ' ...
               '(start %d, remaining %d samples).'], iFrStart, maxSamplesAvail);
    end

    if maxBlocksAvail < iNewNoBlocks
        warning(['AnalyzeRxSig_app: baseband length too short, expected %d blocks, ' ...
                 'only %d available.'], iNewNoBlocks, maxBlocksAvail);
        iNewNoBlocks = maxBlocksAvail;
    end

    totalNeeded = iNewNoBlocks * iNb;

    mFrameRxTp = zeros(totalNeeded, iNoRxAnt);
    mFrameRx   = zeros(iNb, iNewNoBlocks, iNoRxAnt);

    for iC = 1:iNoRxAnt
        seg = mDataRxDem(iC, iFrStart : iFrStart + totalNeeded - 1);
        mFrameRxTp(:,iC) = seg.';
        mFrameRx(:,:,iC) = reshape(seg, iNb, []);
    end

    % Remove CP: [iNfft x iNewNoBlocks x Nr]
    mFrameRxNoCP = mFrameRx;
    mFrameRxNoCP(1:iNg,:,:) = [];

    %% 5. Dispatch to mode-specific RX
    switch lower(mimoModeStr)
        case 'alamouti'
            [mEmpfDataBits, dataMode] = AlamoutiRx_app( ...
                mFrameRxNoCP, params, kanal, len_cInfoBits);

        case {'spatail multiplexing','v-blast','vblast'}
            [mEmpfDataBits, dataMode] = SM_VBLAST_Rx_app( ...
                mFrameRxNoCP, params, kanal, len_cInfoBits);

        case 'eigenmode'
            [mEmpfDataBits, dataMode] = EigenModeRx_app( ...
                mFrameRxNoCP, params, kanal, len_cInfoBits);

        otherwise
            error('AnalyzeRxSig_app: unsupported MIMO mode: %s', mimoModeStr);
    end

    %% 6. Fill common debug info
    data = dataMode;
    data.mRecFrame    = mRecFrame;
    data.fDACFreq     = fDACFreq;
    data.fBBFreq      = fBBFreq;
    data.fCarrFreq    = fCarrFreq;
    data.iNfft        = iNfft;
    data.iModOrd      = iModOrd;
    data.ta_TP        = 1 / fBBFreq;
    data.MetricData   = MetricData;
    data.iNewNoBlocks = iNewNoBlocks;
end

%% ========================================================================
%% Subfunction: SM / V-BLAST RX (preamble extraction + channel est + EQ)
%% ========================================================================
function [mEmpfDataBits, data] = SM_VBLAST_Rx_app(mFrameRxNoCP, params, kanal, len_cInfoBits)
% mFrameRxNoCP : [iNfft x iNewNoBlocks x Nr]

    iNoBlocks     = params.iNoBlocks;
    iNfft         = params.iNfft;
    iNg           = params.iNg;    %#ok<NASGU>
    iNb           = params.iNb;    %#ok<NASGU>
    iModOrd       = params.iModOrd;
    iNoTxAnt      = params.iNoTxAnt;
    iNoRxAnt      = params.iNoRxAnt;
    iNewNoBlocks  = params.iNewNoBlocks;
    iNoSubBlocks  = params.iNoSubBlocks;

    mimoModeStr   = string(params.mimoMode);

    % Subframe structure
    AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);

    % Chu preamble (same as TX)
    vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    vPreambleFreq = 1/sqrt(iNfft) * fft(vPreambleTime, iNfft);

    % 1) Extract preambles
    mPreambleRx = zeros(iNfft, AnzSubFrames, iNoRxAnt, iNoTxAnt);
    for iCTxAnt = 1:iNoTxAnt
        mPreambleRx(:,:,:,iCTxAnt) = ...
            mFrameRxNoCP(1:iNfft, ...
                         iCTxAnt : iNoSubBlocks+iNoTxAnt : iNoSubBlocks*(AnzSubFrames-1)+AnzSubFrames*iNoTxAnt, ...
                         :);
    end

    % 2) Extract data blocks
    v = 0;
    mDataRx = [];
    for ii = 1:AnzSubFrames
        if ii == AnzSubFrames
            mDataRx(:, iNoSubBlocks*(ii-1)+1 : iNoSubBlocks*ii + rem(iNoBlocks,iNoSubBlocks), :) = ...
                mFrameRxNoCP(1:iNfft, ...
                             iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v+rem(iNoBlocks,iNoSubBlocks), ...
                             :);
        else
            mDataRx(:, iNoSubBlocks*(ii-1)+1 : iNoSubBlocks*ii, :) = ...
                mFrameRxNoCP(1:iNfft, ...
                             iNoTxAnt+1+v : iNoTxAnt+iNoSubBlocks+v, ...
                             :);
        end
        v = v + iNoTxAnt + iNoSubBlocks;
    end

    % 3) Tail zero blocks for noise estimation
    if size(mFrameRxNoCP,2) >= 10
        mZeroRx = mFrameRxNoCP(:, end-10+1:end, :);
    else
        mZeroRx = mFrameRxNoCP;
    end

    %% SNR and channel estimation
    mPreambleRxFreq = 1/sqrt(iNfft) * fft(mPreambleRx, iNfft, 1);

    mPowNoise    = zeros(iNfft, iNoRxAnt);
    mPowSigNoisy = zeros(iNfft, iNoRxAnt, iNoTxAnt);
    mSNR         = zeros(iNfft, iNoRxAnt, iNoTxAnt);

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            for k = 1:iNfft
                mPowNoise(k,iCRxAnt) = mean(abs(mZeroRx(k,:,iCRxAnt)).^2);
            end
            mPowSigNoisy(:,iCRxAnt,iCTxAnt) = abs(mPreambleRxFreq(:,end,iCRxAnt,iCTxAnt)).^2;
            mSNR(:,iCRxAnt,iCTxAnt) = ...
                (mPowSigNoisy(:,iCRxAnt,iCTxAnt) - mPowNoise(:,iCRxAnt)) ./ ...
                 max(mPowNoise(:,iCRxAnt), eps);
        end
    end

    % Numerical protection
    badMask = ~isfinite(mSNR) | (mSNR <= 0);
    mSNR(badMask) = 1e6;

    vSNRperSCperBl = zeros(iNfft,1);
    for k = 1:iNfft
        vSNRperSCperBl(k) = mean(mean(squeeze(mSNR(k,:,:))));
    end

    iSNR = mean(vSNRperSCperBl(:));
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    % Channel estimation
    vPreambleFreqCol = vPreambleFreq(:);

    if strcmpi(kanal.schaetzer,'Zero Forcing')
        % LS (ZF)
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                for iCB = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);
                    mCTF(:,iCB,iCRxAnt,iCTxAnt) = rxCol ./ vPreambleFreqCol;
                    mCIR(:,iCB,iCRxAnt,iCTxAnt) = ...
                        sqrt(iNfft) * ifft(mCTF(:,iCB,iCRxAnt,iCTxAnt));
                end
            end
        end
    elseif strcmpi(kanal.schaetzer,'MMSE')
        % MMSE
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                snrCol = mSNR(:,iCRxAnt,iCTxAnt);
                den   = abs(vPreambleFreqCol).^2 + 1./snrCol;
                w     = conj(vPreambleFreqCol) ./ den;
                for iCB = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);
                    mCTF(:,iCB,iCRxAnt,iCTxAnt) = w .* rxCol;
                    mCIR(:,iCB,iCRxAnt,iCTxAnt) = ...
                        sqrt(iNfft) * ifft(mCTF(:,iCB,iCRxAnt,iCTxAnt));
                end
            end
        end
    else
        error('SM_VBLAST_Rx_app: unknown channel estimator: %s', kanal.schaetzer);
    end

    %% MIMO equalization
    % Data freq: [SC x Block x Rx] -> [Rx x SC x Block]
    mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);
    mDataRxFreq = permute(mDataRxFreq, [3,1,2]);
    mDataRxEq   = zeros(iNoTxAnt, iNfft, iNoBlocks);

    if strcmpi(mimoModeStr,'v-blast') || strcmpi(mimoModeStr,'vblast')
        % V-BLAST detector
        mDataRxEq = VBlastRx_app( ...
            mDataRx, ...
            mCTF, ...
            struct( ...
                'iNoBlocks',    iNoBlocks, ...
                'iNfft',        iNfft, ...
                'iNoTxAnt',     iNoTxAnt, ...
                'iNoRxAnt',     iNoRxAnt, ...
                'iNoSubBlocks', iNoSubBlocks ...
            ), ...
            iSNR, ...
            kanal ...
        );
    else
        % Linear ZF / MMSE for spatial multiplexing
        if strcmpi(kanal.entzerrer,'Zero Forcing')
            for iSF = 1:AnzSubFrames
                if iSF == AnzSubFrames
                    iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
                else
                    iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
                end
                for iCB = iCB_range
                    for k = 1:iNfft
                        mH = squeeze(mCTF(k,iSF,:,:));   % [Rx x Tx]
                        mDataRxEq(:,k,iCB) = sqrt(iNoTxAnt) * pinv(mH) * mDataRxFreq(:,k,iCB);
                    end
                end
            end

        elseif strcmpi(kanal.entzerrer,'MMSE')
            for iSF = 1:AnzSubFrames
                if iSF == AnzSubFrames
                    iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks + rem(iNoBlocks,iNoSubBlocks);
                else
                    iCB_range = (iSF-1)*iNoSubBlocks+1 : iSF*iNoSubBlocks;
                end
                for iCB = iCB_range
                    for k = 1:iNfft
                        mH = squeeze(mCTF(k,iSF,:,:));
                        A  = mH' * mH + eye(iNoTxAnt) * (1./iSNR);
                        if any(~isfinite(A(:))) || rcond(A) < 1e-10
                            A = mH' * mH + eye(iNoTxAnt) * 1e-3;
                        end
                        mDataRxEq(:,k,iCB) = ...
                            sqrt(iNoTxAnt) * ( A \ (mH' * mDataRxFreq(:,k,iCB)) );
                    end
                end
            end
        else
            error('SM_VBLAST_Rx_app: unknown equalizer: %s', kanal.entzerrer);
        end
    end

    % Maximal ratio combining for single TX antenna (kept from legacy code)
    if iNoTxAnt == 1
        mDataRxEqTp    = mDataRxEq;
        mWeightFactors = zeros(iNfft, iNoBlocks);
        for iCB = 1:AnzSubFrames
            for k = 1:iNfft
                mWeightFactors(k,iCB) = norm(squeeze(mCTF(k,iCB,:,:)));
            end
            mWeightFactors(:,iCB) = mWeightFactors(:,iCB) / max(sum(mWeightFactors(:,iCB)),eps);
        end
        for iCB = 1:iNoBlocks
            for k = 1:iNfft
                mDataRxEqTp(:,k,iCB) = mDataRxEq(:,k,iCB) * mWeightFactors(k,iCB);
            end
        end
        mDataRxEq = mDataRxEqTp;
    end

    %% Symbol -> bits
    for iCTxAnt = 1:iNoTxAnt
        mDataRxDet(:,:,iCTxAnt) = squeeze(mDataRxEq(iCTxAnt,:,:));
    end
    vDataRxDet = reshape(mDataRxDet, 1, iNfft*iNoBlocks*iNoTxAnt);

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

    vBitsRxDetTp  = [vInfoRxDetTp(:); vBitsRxDet];
    mEmpfDataBits = vBitsRxDetTp.';

    %% Pack debug data
    data             = struct();
    if exist('mCTF','var')
        data.kanalUeb = mCTF;
    else
        data.kanalUeb = [];
    end
    if exist('mCIR','var')
        data.kanalimpuls = mCIR;
    else
        data.kanalimpuls = [];
    end
    data.mDataRxEq    = mDataRxEq;
    data.AnzSubFrames = AnzSubFrames;
end