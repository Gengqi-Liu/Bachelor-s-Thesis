function [mFrameTxCar, meta] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime)
%GENERATETX_SM_APP
% Tx generation for Spatial Multiplexing / V-BLAST:
%   - Per-antenna OFDM
%   - Frame: preambles + data subframes
%   - CP, sync, resampling and upconversion

    % Unpack parameters
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iNg          = params.iNg;
    iNb          = params.iNb;
    iNoTxAnt     = params.iNoTxAnt;
    iNoSubBlocks = params.iNoSubBlocks;

    fBBFreq      = params.fBBFreq;
    fDACFreq     = params.fDACFreq;
    fCarrFreq    = params.fCarrFreq;

    meta = struct();

    % Per-antenna power normalization
    mDataTxFreq = 1/sqrt(iNoTxAnt) * mDataTxFreq;

    % OFDM modulation (frequency -> time)
    mDataTxTime = sqrt(iNfft) * ifft(mDataTxFreq, iNfft, 1);   % [iNfft x iNoBlocks x iNoTxAnt]

    % ---- Build frame: Preamble + subframes of data ----
    AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);
    mFrame       = zeros(iNfft, iNoBlocks + AnzSubFrames*iNoTxAnt, iNoTxAnt);

    m = 0;
    n = 0;
    b = iNoSubBlocks;

    for ii = 1:(2*AnzSubFrames)
        if rem(ii,2) ~= 0
            % Odd step: insert preamble (one per Tx antenna)
            ind = (1+m):(iNoTxAnt+m);
            for iCTxAnt = 1:iNoTxAnt
                mFrame(:,ind(iCTxAnt),iCTxAnt) = vPreambleTime;
            end
            m = m + b + iNoTxAnt;

        elseif ii == 2*AnzSubFrames
            % Last segment: may contain tail blocks
            ind = (iNoTxAnt+1+n):(iNoTxAnt + b + n + rem(iNoBlocks,b));
            mFrame(:,ind,:) = mDataTxTime(:, b*(ii/2-1)+1 : b*ii/2 + rem(iNoBlocks,b), :);
        else
            % Even step: full subframe of data blocks
            ind = (iNoTxAnt+1+n):(iNoTxAnt + b + n);
            mFrame(:,ind,:) = mDataTxTime(:, b*(ii/2-1)+1 : b*ii/2, :);
            n = n + b + iNoTxAnt;
        end
    end

    % A) Append empty OFDM blocks at tail
    mFrameTxTp = cat(2, mFrame, zeros(iNfft,15 , iNoTxAnt));

    % B) Add cyclic prefix (baseband) - robust for iNg > iNfft
    L = size(mFrameTxTp,1);   % useful symbol length (should be iNfft)
    
    if iNg <= 0
        cpPart = zeros(0, size(mFrameTxTp,2), size(mFrameTxTp,3)); % empty CP
    elseif iNg <= L
        cpPart = mFrameTxTp(end-iNg+1:end,:,:);
    else
        % iNg > L: build CP by cyclically repeating the OFDM symbol tail
        % idx has length iNg and values always in 1..L
        idx = mod((L - iNg):(L - 1), L) + 1;
        cpPart = mFrameTxTp(idx,:,:);
    end
    
    mFrameBB = [cpPart; mFrameTxTp];   % [(iNg+iNfft) x iNewNoBlocks x iNoTxAnt]

    szBB              = size(mFrameBB);
    iNewNoBlocks      = szBB(2);
    meta.iNewNoBlocks = iNewNoBlocks;

    % Serialize baseband payload
    payloadBB = reshape(mFrameBB, iNb * iNewNoBlocks, iNoTxAnt);  % [Ns_payload x iNoTxAnt]
    meta.payloadLen = size(payloadBB,1);

    % C) Insert Schmidl & Cox sync block in baseband
    if exist('buildSchmidlCoxSync_app','file')
        [syncNoCp, Lsync] = buildSchmidlCoxSync_app(iNfft);
    else
        error('generateTx_SM_app: buildSchmidlCoxSync_app.m not found on path.');
    end
    syncBlock = repmat(syncNoCp, 1, iNoTxAnt);    % [iNfft x iNoTxAnt]

    headZeros = zeros(1000, iNoTxAnt);
    %tailZeros = zeros(1, iNoTxAnt);

    % Final baseband: silence + sync + payload + silence
    mFrameTxBB = [headZeros; syncBlock; payloadBB];

    meta.Lsync        = Lsync;
    meta.syncLen      = length(syncNoCp);
    meta.headLen      = size(headZeros,1);
    %meta.tailLen      = size(tailZeros,1);
    meta.basebandLen  = size(mFrameTxBB,1);

    % D) Resample to DAC rate
    mFrameTxUp = zeros( ceil(size(mFrameTxBB,1) * fDACFreq / fBBFreq), iNoTxAnt );
    for iCAnt = 1:iNoTxAnt
        mFrameTxUp(:,iCAnt) = resample(mFrameTxBB(:,iCAnt), fDACFreq, fBBFreq);
    end

    % E) Upconvert to carrier
    Ns = size(mFrameTxUp,1);
    t  = (0:Ns-1).' / fDACFreq;
    mFrameTxCar = zeros(Ns, iNoTxAnt);
    for iCAnt = 1:iNoTxAnt
        mFrameTxCar(:,iCAnt) = real(mFrameTxUp(:,iCAnt) .* exp(1j*2*pi*fCarrFreq*t));
    end

    % F) Amplitude normalization
    maxVal = max(abs(mFrameTxCar), [], 'all');
    if maxVal > 0
        mFrameTxCar = 1.5/maxVal * mFrameTxCar;
    end

    % G) Estimated playback time
    wavTime_sec = size(mFrameTxCar,1) / fDACFreq;
    wavTime_ms  = wavTime_sec * 1000 * 1.3;
    wavTime     = ceil(wavTime_ms)/1000;
    meta.wavTime = wavTime;
end