function [U, S, V, chanMeta] = EigenMode_ChanEst_app(rxPreamble, params, metaPre)
% EigenMode_ChanEst_app
%
% Estimate frequency-domain MIMO channel from EigenMode preamble and compute
% per-subcarrier SVD of (H^H H). V is used for Tx precoding.
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

    % ---- Legacy sync (Chu) via EstFrameStart ----
    PglobalTmp = struct( ...
        'iNoBlocks', iNoBlocks, ...
        'iNfft',     iNfft, ...
        'iNg',       iNg, ...
        'iNb',       iNb, ...
        'iNoTxAnt',  iNoTxAnt, ...
        'iNoRxAnt',  iNoRxAnt, ...
        'fBBFreq',   fBBFreq, ...
        'fDACFreq',  fDACFreq, ...
        'fCarrFreq', fCarrFreq );

    vFrStart   = [];
    vCarOff    = [];
    MetricData = cell(1, iNoRxAnt);

    for iC = 1:iNoRxAnt
        [vFrStartTp, fCarOffTp, MetricData{iC}] = EstFrameStart(mDataRxDem(iC,:), PglobalTmp);
        vFrStart = [vFrStart, vFrStartTp]; %#ok<AGROW>
        vCarOff  = [vCarOff,  fCarOffTp];  %#ok<AGROW>
    end

    iFrStart = max(1, min(vFrStart) - 4);

    vCarOff = vCarOff(:).';
    fCarOff = mean(vCarOff);

    % ---- Per-Rx coarse CFO correction ----
    Ns = size(mDataRxDem,2);
    n  = 1:Ns;
    for iC = 1:iNoRxAnt
        mDataRxDem(iC,:) = mDataRxDem(iC,:) .* exp(-1j*2*pi/iNfft * vCarOff(iC) * n);
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

    mPreambleRx(1:iNg,:,:) = [];  % -> [iNfft x iNewNoBlocksPre x Nr]

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

    mPreambleRxFreq1 = reshape(mPreambleRxFreq, iNfft, iNoTxAnt, N, iNoRxAnt);

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
    chanMeta.iFrStart        = iFrStart;
    chanMeta.fCarOff         = fCarOff;
end