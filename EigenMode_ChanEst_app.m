function [U, S, V, chanMeta] = EigenMode_ChanEst_app(rxPreamble, params, metaPre)
% EigenMode_ChanEst_app
%
% Estimate frequency-domain MIMO channel from EigenMode preamble and compute
% per-subcarrier SVD of mean(H^H H). V is used for Tx precoding.
%
% Input:
%   rxPreamble : [Nr x Nsamp] recorded HF signal (transpose of audio I/O)
%   params     : struct with fields:
%       iNoBlocks, iNfft, iNg, iNb, iNoTxAnt, iNoRxAnt, fBBFreq, fDACFreq, fCarrFreq
%   metaPre    : struct with fields:
%       iNewNoBlocksPre (= N * iNoTxAnt), N
%
% Output:
%   U,S,V      : [Nt x Nt x iNfft]
%   chanMeta   : struct with debug fields (mCTF, mCIR, iFrStart, fCarOff, MetricData, ...)

    % ---- Unpack params ----
    iNoBlocks       = params.iNoBlocks;
    iNfft           = params.iNfft;
    iNg             = params.iNg;
    iNb             = params.iNb;
    iNoTxAnt        = params.iNoTxAnt;
    iNoRxAnt        = params.iNoRxAnt;

    fBBFreq         = params.fBBFreq;
    fDACFreq        = params.fDACFreq;
    fCarrFreq       = params.fCarrFreq;

    iNewNoBlocksPre = metaPre.iNewNoBlocksPre;
    N               = metaPre.N;

    chanMeta = struct();

    % rxPreamble: [Nr x Nsamp]
    mRecPreamble = rxPreamble;

    % ---- Downconvert to baseband ----
    vBB1 = DeModulateSignal_app(mRecPreamble(1,:), fBBFreq, fDACFreq, fCarrFreq);
    Lbb  = numel(vBB1);

    mDataRxDem      = zeros(iNoRxAnt, Lbb);
    mDataRxDem(1,:) = vBB1;

    for iC = 2:iNoRxAnt
        vBB = DeModulateSignal_app(mRecPreamble(iC,:), fBBFreq, fDACFreq, fCarrFreq);
        if numel(vBB) ~= Lbb
            L = min(Lbb, numel(vBB));
            mDataRxDem(iC,1:L) = vBB(1:L);
        else
            mDataRxDem(iC,:) = vBB;
        end
    end

    % ---- Schmidl & Cox timing + coarse CFO (uses EstFrameStart_app) ----
    if ~exist('EstFrameStart_app','file')
        error('EigenMode_ChanEst_app: EstFrameStart_app.m not found on path.');
    end

    vFrStart   = zeros(1, iNoRxAnt);
    vCarOff    = zeros(1, iNoRxAnt);
    MetricData = cell(1, iNoRxAnt);

    for iC = 1:iNoRxAnt
        [vFrStart(iC), vCarOff(iC), MetricData{iC}] = EstFrameStart_app( ...
            mDataRxDem(iC,:), iNfft, iNg);
    end

    % Use earliest sync point across Rx channels
    iFrSync  = min(vFrStart);

    % Data/preamble blocks start AFTER the S&C sync symbol (length iNfft, no CP)
    iFrStart = iFrSync + iNfft;
    iFrStart = max(1, min(iFrStart, Lbb));

    % Coarse CFO compensation (normalized cycles/sample)
    cfoVals = vCarOff(~isnan(vCarOff) & isfinite(vCarOff));
    if ~isempty(cfoVals)
        fCfo = median(cfoVals);  % robust across channels
    else
        fCfo = 0;
    end

    if isfinite(fCfo) && abs(fCfo) > 0
        n = 0:(Lbb-1);
        rot = exp(-1j*2*pi*fCfo*n);
        for iC = 1:iNoRxAnt
            mDataRxDem(iC,:) = mDataRxDem(iC,:) .* rot;
        end
    end

    % ---- Slice preamble blocks and remove CP ----
    totalNeeded = iNewNoBlocksPre * iNb;
    maxAvail    = Lbb - iFrStart + 1;

    if maxAvail < totalNeeded
        error('EigenMode_ChanEst_app: insufficient preamble length (need %d, have %d).', ...
              totalNeeded, maxAvail);
    end

    mPreambleRxTp = zeros(totalNeeded, iNoRxAnt);
    mPreambleRx   = zeros(iNb, iNewNoBlocksPre, iNoRxAnt);

    for iC = 1:iNoRxAnt
        idx = iFrStart + (0:totalNeeded-1);
        mPreambleRxTp(:,iC) = mDataRxDem(iC, idx).';
        mPreambleRx(:,:,iC) = reshape(mPreambleRxTp(:,iC), iNb, iNewNoBlocksPre);
    end

    % Remove CP: -> [iNfft x iNewNoBlocksPre x Nr]
    mPreambleRx(1:iNg,:,:) = [];

    % ---- ZF channel estimation in frequency domain ----
    if exist('ChuSeq','file')
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft);
    else
        warning('EigenMode_ChanEst_app: ChuSeq not found, using random preamble.');
        vPreambleTime = sqrt(iNoTxAnt) * (randn(1,iNfft) + 1j*randn(1,iNfft))/sqrt(2);
    end
    vPreambleTime = vPreambleTime(:).';
    vPreambleFreq = 1/sqrt(iNfft) * fft(vPreambleTime, iNfft);

    mPreambleRxFreq = 1/sqrt(iNfft) * fft(mPreambleRx, iNfft, 1);

    if iNewNoBlocksPre ~= iNoTxAnt * N
        error('EigenMode_ChanEst_app: iNewNoBlocksPre (%d) must equal iNoTxAnt*N (%d).', ...
              iNewNoBlocksPre, iNoTxAnt*N);
    end

    % Reshape: [SC x (Nt*N) x Nr] -> [SC x Nt x N x Nr]
    mPreambleRxFreq1 = reshape(mPreambleRxFreq, iNfft, iNoTxAnt, N, iNoRxAnt);

    % mCTF: [SC x Nt x Nr x N]
    mCTF = zeros(iNfft, iNoTxAnt, iNoRxAnt, N);
    mCIR = zeros(iNfft, iNoTxAnt, iNoRxAnt, N);

    vPreambleFreq_col = vPreambleFreq(:);

    for iCTxAnt = 1:iNoTxAnt
        for iCRxAnt = 1:iNoRxAnt
            for iRep = 1:N
                rxCol = mPreambleRxFreq1(:, iCTxAnt, iRep, iCRxAnt);
                Hcol  = rxCol ./ vPreambleFreq_col;
                mCTF(:,iCTxAnt,iCRxAnt,iRep) = Hcol;
                mCIR(:,iCTxAnt,iCRxAnt,iRep) = sqrt(iNfft) * ifft(Hcol);
            end
        end
    end

    % ---- Per-subcarrier SVD of mean(H^H H) ----
    U = zeros(iNoTxAnt, iNoTxAnt, iNfft);
    S = zeros(iNoTxAnt, iNoTxAnt, iNfft);
    V = zeros(iNoTxAnt, iNoTxAnt, iNfft);

    for iSc = 1:iNfft
        H2_mean = zeros(iNoTxAnt, iNoTxAnt);

        for iRep = 1:N
            H_temp = zeros(iNoRxAnt, iNoTxAnt);
            for iCTxAnt = 1:iNoTxAnt
                H_temp(:,iCTxAnt) = squeeze(mCTF(iSc,iCTxAnt,:,iRep));
            end
            H2_mean = H2_mean + (H_temp' * H_temp);
        end

        if N > 1
            H2_mean = H2_mean / N;
        end

        [u,s,v] = svd(H2_mean);
        U(:,:,iSc) = u;
        S(:,:,iSc) = s;
        V(:,:,iSc) = v;
    end

    % ---- Debug output ----
    chanMeta.mCTF            = mCTF;
    chanMeta.mCIR            = mCIR;
    chanMeta.mPreambleRx     = mPreambleRx;
    chanMeta.mPreambleRxFreq = mPreambleRxFreq;
    chanMeta.MetricData      = MetricData;
    chanMeta.vFrStart        = vFrStart;
    chanMeta.vCarOff         = vCarOff;
    chanMeta.iFrStart        = iFrStart;
    chanMeta.fCarOff         = fCfo;
end
