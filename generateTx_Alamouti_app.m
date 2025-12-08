function [mFrameTxCar, meta] = generateTx_Alamouti_app(params, mDataTxFreq, vPreambleTime)
%GENERATETX_ALAMOUTI_APP
% Tx generation for Alamouti 2x?:
%   - Frequency-domain Alamouti mapping (2 time slots per OFDM block)
%   - Frame structure with preambles and nulls
%   - CP, sync block, resampling and upconversion

    % --- Unpack parameters ---
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iNg          = params.iNg;
    iNb          = params.iNb;
    iNoTxAnt     = params.iNoTxAnt;
    fBBFreq      = params.fBBFreq;
    fDACFreq     = params.fDACFreq;
    fCarrFreq    = params.fCarrFreq;

    if iNoTxAnt ~= 2
        error('generateTx_Alamouti_app: currently only 2-Tx Alamouti is implemented (iNoTxAnt must be 2).');
    end

    meta = struct();

    %% 1) Frequency-domain Alamouti encoding (2x time slots per block)
    % Original mapping:
    %   Time slot 1: [s1, s2]
    %   Time slot 2: [-conj(s2), conj(s1)]

    nBlocksAlam = 2 * iNoBlocks;    % logical OFDM blocks after Alamouti

    MimoDataTxFreq = zeros(iNfft, nBlocksAlam, iNoTxAnt);

    for iCB = 1:iNoBlocks
        % Slot 1: direct symbols
        MimoDataTxFreq(:, 2*iCB-1, :) = mDataTxFreq(:, iCB, :);

        % Slot 2: [-s2* , s1*]
        MimoDataTxFreq(:, 2*iCB, 1) = conj(-mDataTxFreq(:, iCB, 2));
        MimoDataTxFreq(:, 2*iCB, 2) = conj( mDataTxFreq(:, iCB, 1));
    end

    % Power normalization over Tx antennas
    MimoDataTxFreq = 1/sqrt(iNoTxAnt) * MimoDataTxFreq;

    % Frequency -> time
    MimoDataTxTime = sqrt(iNfft) * ifft(MimoDataTxFreq, iNfft, 1);  % [iNfft x 2*iNoBlocks x 2]

    %% 2) Frame structure with preambles and null columns
    % Layout per logical block (Nt = 2, nColsPerBlock = 5):
    %   columns: [Tx1 preamble, slot1, slot2, Tx2 preamble, null]

    nColsPerBlock = iNoTxAnt + 3;   % = 5
    nFrameCols    = nColsPerBlock * iNoBlocks;

    mFrameTxTp = zeros(iNfft, nFrameCols, iNoTxAnt);

    for iCB = 1:iNoBlocks
        vInd = (nColsPerBlock*(iCB-1)+1) : (nColsPerBlock*iCB);  % length 5

        % Insert preambles
        % Nt=2:
        %   Tx1 preamble at vInd(1) + 0
        %   Tx2 preamble at vInd(2) + 2
        for iCAnt = 1:iNoTxAnt
            colP = vInd(iCAnt) + (iCAnt-1)*2;
            mFrameTxTp(:, colP, iCAnt) = vPreambleTime;
        end

        % Insert two Alamouti time slots in the middle
        mFrameTxTp(:, vInd(end)-3 : vInd(end)-2, :) = ...
            MimoDataTxTime(:, 2*iCB-1 : 2*iCB, :);

        % Last column vInd(end) remains zero (null)
    end

    % Append additional empty OFDM blocks (tail)
    mFrameTxTp = cat(2, mFrameTxTp, zeros(iNfft, 20, iNoTxAnt));

    %% 3) Add cyclic prefix
    mFrameBB = [mFrameTxTp(end-iNg+1:end,:,:); mFrameTxTp];  % [iNb x iNewNoBlocks x 2]

    szBB              = size(mFrameBB);
    iNewNoBlocks      = szBB(2);
    meta.iNewNoBlocks = iNewNoBlocks;

    % Serialize baseband payload
    payloadBB        = reshape(mFrameBB, iNb * iNewNoBlocks, iNoTxAnt);
    meta.payloadLen  = size(payloadBB,1);

    %% 4) Add Schmidl & Cox sync block and head/tail zeros
    if exist('buildSchmidlCoxSync_app','file')
        [syncNoCp, Lsync] = buildSchmidlCoxSync_app(iNfft);
    else
        error('generateTx_Alamouti_app: buildSchmidlCoxSync_app.m not found on path.');
    end

    syncBlock = repmat(syncNoCp, 1, iNoTxAnt);   % [iNfft x 2]

    headZeros = zeros(1000, iNoTxAnt);
    tailZeros = zeros(3000, iNoTxAnt);

    mFrameTxBB = [headZeros; syncBlock; payloadBB; tailZeros];

    meta.Lsync       = Lsync;
    meta.syncLen     = length(syncNoCp);
    meta.headLen     = size(headZeros,1);
    meta.tailLen     = size(tailZeros,1);
    meta.basebandLen = size(mFrameTxBB,1);

    %% 5) Resample to DAC rate
    nBB  = size(mFrameTxBB,1);
    nUp  = ceil(nBB * fDACFreq / fBBFreq);
    mFrameTxUp = zeros(nUp, iNoTxAnt);

    for iCAnt = 1:iNoTxAnt
        mFrameTxUp(:,iCAnt) = resample(mFrameTxBB(:,iCAnt), fDACFreq, fBBFreq);
    end

    %% 6) Upconvert to carrier frequency
    Ns = size(mFrameTxUp,1);
    t  = (0:Ns-1).' / fDACFreq;

    mFrameTxCar = zeros(Ns, iNoTxAnt);
    for iCAnt = 1:iNoTxAnt
        mFrameTxCar(:,iCAnt) = real(mFrameTxUp(:,iCAnt) .* exp(1j*2*pi*fCarrFreq*t));
    end

    %% 7) Amplitude normalization and duration estimate
    maxVal = max(abs(mFrameTxCar), [], 'all');
    if maxVal > 0
        mFrameTxCar = 1.5 / maxVal * mFrameTxCar;
    end

    wavTime_sec = size(mFrameTxCar,1) / fDACFreq;
    wavTime_ms  = wavTime_sec * 1000 * 1.3;
    wavTime     = ceil(wavTime_ms) / 1000;
    meta.wavTime = wavTime;
end