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
    iNfft         = params.iNfft;
    iNg           = params.iNg;
    iNb           = params.iNb;
    iModOrd       = params.iModOrd;
    iNoRxAnt      = params.iNoRxAnt;
    iNewNoBlocks  = params.iNewNoBlocks;

    fBBFreq       = params.fBBFreq;
    fDACFreq      = params.fDACFreq;
    fCarrFreq     = params.fCarrFreq;

    mimoModeStr   = string(params.mimoMode);

    if isfield(kanal,'len_cInfoBits') && ~isempty(kanal.len_cInfoBits)
        len_cInfoBits = kanal.len_cInfoBits;
    else
        len_cInfoBits = 0;
    end

    mEmpfDataBits = [];
    data          = struct();

    %% 1.5 透传数据类型（关键：图像模式需要在RX端分支）
    % 优先 kanal.DatenTyp，其次 params.DatenTyp（如果你有）
    if ~isfield(kanal,'DatenTyp') || isempty(kanal.DatenTyp)
        if isfield(params,'DatenTyp') && ~isempty(params.DatenTyp)
            kanal.DatenTyp = params.DatenTyp;
        end
    end

    % 也把 iModOrd 放进 kanal，给 bits2bild_app/图像header重复解析用
    if ~isfield(kanal,'iModOrd') || isempty(kanal.iModOrd)
        kanal.iModOrd = iModOrd;
    end

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
    iFrStart = iFrSync + iNfft;
    iFrStart = max(1, min(iFrStart, Lbb));

    % --- CFO compensation (fCarOffset is normalized cycles/sample) ---
    if ~isempty(vCarOff)
        cfoVals = vCarOff(~isnan(vCarOff) & isfinite(vCarOff));
        if ~isempty(cfoVals)
            fCfo = median(cfoVals);  % robust across Rx channels
            if isfinite(fCfo) && abs(fCfo) > 0
                n = 0:(Lbb-1);
                rot = exp(-1j*2*pi*fCfo*n);  % cycles/sample -> rad/sample
                for iC = 1:iNoRxAnt
                    mDataRxDem(iC,:) = mDataRxDem(iC,:) .* rot;
                end
            end
        end
    end

    %% 4. Cut OFDM blocks and remove CP
    maxSamplesAvail = Lbb - iFrStart + 1;
    maxBlocksAvail  = floor(maxSamplesAvail / iNb);

    if maxBlocksAvail <= 0
        error(['AnalyzeRxSig_app: not enough samples for one OFDM block ' ...
               '(start %d, remaining %d samples).'], iFrStart, maxSamplesAvail);
    end

    if maxBlocksAvail < iNewNoBlocks
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

        case {'spatial multiplexing','v-blast','vblast'}
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
    iModOrd       = params.iModOrd;
    iNoTxAnt      = params.iNoTxAnt;
    iNoRxAnt      = params.iNoRxAnt;
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
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                for iCB = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);
                    mCTF(:,iCB,iCRxAnt,iCTxAnt) = rxCol ./ vPreambleFreqCol;
                    mCIR(:,iCB,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(mCTF(:,iCB,iCRxAnt,iCTxAnt));
                end
            end
        end
    elseif strcmpi(kanal.schaetzer,'MMSE')
        for iCTxAnt = 1:iNoTxAnt
            for iCRxAnt = 1:iNoRxAnt
                snrCol = mSNR(:,iCRxAnt,iCTxAnt);
                den   = abs(vPreambleFreqCol).^2 + 1./snrCol;
                w     = conj(vPreambleFreqCol) ./ den;
                for iCB = 1:AnzSubFrames
                    rxCol = mPreambleRxFreq(:,iCB,iCRxAnt,iCTxAnt);
                    mCTF(:,iCB,iCRxAnt,iCTxAnt) = w .* rxCol;
                    mCIR(:,iCB,iCRxAnt,iCTxAnt) = sqrt(iNfft) * ifft(mCTF(:,iCB,iCRxAnt,iCTxAnt));
                end
            end
        end
    else
        error('SM_VBLAST_Rx_app: unknown channel estimator: %s', kanal.schaetzer);
    end

    %% MIMO equalization
    mDataRxEq   = zeros(iNoTxAnt, iNfft, iNoBlocks);

    if strcmpi(mimoModeStr,'v-blast') || strcmpi(mimoModeStr,'vblast')
        mDataRxEq = VBlastRx_app( ...
            mDataRx, ...
            mCTF, ...
            struct( ...
                'iNoBlocks',    iNoBlocks, ...
                'iNfft',        iNfft, ...
                'iNoTxAnt',     iNoTxAnt, ...
                'iNoRxAnt',     iNoRxAnt, ...
                'iNoSubBlocks', iNoSubBlocks, ...
                'iModOrd',      iModOrd ...
            ), ...
            iSNR, ...
            kanal ...
        );
    else
        % 线性ZF/MMSE（保留你原实现）
        mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);
        mDataRxFreq = permute(mDataRxFreq, [3,1,2]); % [Rx x SC x Block]

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
                        mDataRxEq(:,k,iCB) = sqrt(iNoTxAnt) * ( A \ (mH' * mDataRxFreq(:,k,iCB)) );
                    end
                end
            end
        else
            error('SM_VBLAST_Rx_app: unknown equalizer: %s', kanal.entzerrer);
        end
    end

    %% Symbol -> bits  （这里是关键改动：图像模式全QAM）
    for iCTxAnt = 1:iNoTxAnt
        mDataRxDet(:,:,iCTxAnt) = squeeze(mDataRxEq(iCTxAnt,:,:));
    end
    vDataRxDet = reshape(mDataRxDet, 1, iNfft*iNoBlocks*iNoTxAnt);
    vDataRxDet(~isfinite(vDataRxDet)) = 0;

    isImage = false;
    if isfield(kanal,'DatenTyp') && ~isempty(kanal.DatenTyp)
        dt = lower(string(kanal.DatenTyp));
        isImage = (dt=="bild" || dt=="image");
    end

    if isImage
        % 图像：整段符号统一按 2^iModOrd 解调成 bit（header也在里面）
        symIdx  = qamdemod(vDataRxDet, 2^iModOrd, 'UnitAveragePower', true);
        bitsMat = de2bi(symIdx, iModOrd, 'left-msb').';
        mEmpfDataBits = bitsMat(:).';   % row bits
    else
        % 文本：保留原来的"BPSK header + QAM payload"
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
    end

    %% Pack debug data
    data = struct();
    data.kanalUeb     = mCTF;
    data.kanalimpuls  = mCIR;
    data.mDataRxEq    = mDataRxEq;
    data.AnzSubFrames = AnzSubFrames;
end
