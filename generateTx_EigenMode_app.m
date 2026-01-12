function [mFrameTxCar, meta] = generateTx_EigenMode_app(params, mDataTxFreq, vPreambleTime)
%GENERATETX_EIGENMODE_APP
% Tx generation for EigenMode:
%   - Frequency-domain precoding using V1(:,:,k)
%   - Standard OFDM framing, CP, sync, upsampling and upconversion

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

    % ---- Load / construct V1 (precoding matrices) ----
    if ~isfield(params, 'V1') || isempty(params.V1)
        warning(['generateTx_EigenMode_app: params.V1 not provided. ', ...
                 'Using identity precoder (EigenMode degenerates to SM).']);
        V1 = repmat(eye(iNoTxAnt), 1, 1, iNfft);   % [Nt x Nt x Nfft]
    else
        V1 = params.V1;
        szV = size(V1);
        if numel(szV) ~= 3 || szV(1) ~= iNoTxAnt || szV(2) ~= iNoTxAnt || szV(3) ~= iNfft
            error('params.V1 has incompatible size. Expected [iNoTxAnt x iNoTxAnt x iNfft].');
        end
    end
    
    meta = struct();
    meta.V1 = V1;

    % ==== 1. Frequency-domain EigenMode precoding ====
    % mDataTxFreq: [iNfft x iNoBlocks x iNoTxAnt]
    mDataTxFreqEq = zeros(size(mDataTxFreq));
    for k = 1:iNoBlocks
        for iSc = 1:iNfft
            sVec = squeeze(mDataTxFreq(iSc, k, :));          % [iNoTxAnt x 1]
            mDataTxFreqEq(iSc, k, :) = V1(:,:,iSc) * sVec;   % [iNoTxAnt x 1]
        end
    end

    mDataTxFreq = mDataTxFreqEq;

    % ==== 2. OFDM processing (same as SM) ====

    % Power normalization over Tx antennas
    mDataTxFreq = 1/sqrt(iNoTxAnt) * mDataTxFreq;

    % IFFT: frequency -> time
    mDataTxTime = sqrt(iNfft) * ifft(mDataTxFreq, iNfft, 1);

    % ---- Insert preambles and build OFDM frame ----
    AnzSubFrames = floor(iNoBlocks / iNoSubBlocks);
    mFrame       = zeros(iNfft, iNoBlocks + AnzSubFrames*iNoTxAnt, iNoTxAnt);

    m = 0;
    n = 0;
    b = iNoSubBlocks;

    for ii = 1:(2*AnzSubFrames)
        if rem(ii,2) ~= 0
            % Odd step: insert preamble for each Tx antenna
            ind = (1+m):(iNoTxAnt+m);
            for iCTxAnt = 1:iNoTxAnt
                mFrame(:,ind(iCTxAnt),iCTxAnt) = vPreambleTime;
            end
            m = m + b + iNoTxAnt;

        elseif ii == 2*AnzSubFrames
            % Last step: possibly partial tail
            ind = (iNoTxAnt+1+n):(iNoTxAnt + b + n + rem(iNoBlocks,b));
            mFrame(:,ind,:) = mDataTxTime(:, b*(ii/2-1)+1 : b*ii/2 + rem(iNoBlocks,b), :);
        else
            % Even step: full sub-blocks
            ind = (iNoTxAnt+1+n):(iNoTxAnt + b + n);
            mFrame(:,ind,:) = mDataTxTime(:, b*(ii/2-1)+1 : b*ii/2, :);
            n = n + b + iNoTxAnt;
        end
    end

    % ---- A. Append tail of empty OFDM symbols ----
    mFrameTxTp = cat(2, mFrame, zeros(iNfft, 15, iNoTxAnt));

    % ---- B. Add cyclic prefix (baseband) ----
    mFrameBB = [mFrameTxTp(end-iNg+1:end,:,:); mFrameTxTp];  % [iNb x iNewNoBlocks x iNoTxAnt]

    szBB              = size(mFrameBB);
    iNewNoBlocks      = szBB(2);
    meta.iNewNoBlocks = iNewNoBlocks;

    % Serialize baseband payload
    payloadBB = reshape(mFrameBB, iNb * iNewNoBlocks, iNoTxAnt);
    meta.payloadLen = size(payloadBB,1);

    % ---- C. Insert Schmidl & Cox sync block ----
    if exist('buildSchmidlCoxSync_app','file')
        [syncNoCp, Lsync] = buildSchmidlCoxSync_app(iNfft);
    else
        error('buildSchmidlCoxSync_app.m not found on path.');
    end
    syncBlock = repmat(syncNoCp, 1, iNoTxAnt);    % [iNfft x iNoTxAnt]

    headZeros = zeros(1000, iNoTxAnt);
    %tailZeros = zeros(1000, iNoTxAnt);

    mFrameTxBB = [headZeros; syncBlock; payloadBB];

    meta.Lsync        = Lsync;
    meta.syncLen      = length(syncNoCp);
    meta.headLen      = size(headZeros,1);
   % meta.tailLen      = size(tailZeros,1);
    meta.basebandLen  = size(mFrameTxBB,1);

    % ---- D. Resample to DAC rate ----
    mFrameTxUp = zeros( ceil(size(mFrameTxBB,1) * fDACFreq / fBBFreq), iNoTxAnt );
    for iCAnt = 1:iNoTxAnt
        mFrameTxUp(:,iCAnt) = resample(mFrameTxBB(:,iCAnt), fDACFreq, fBBFreq);
    end

    % ---- E. Upconversion to carrier ----
    Ns = size(mFrameTxUp,1);
    t  = (0:Ns-1).' / fDACFreq;
    mFrameTxCar = zeros(Ns, iNoTxAnt);
    for iCAnt = 1:iNoTxAnt
        mFrameTxCar(:,iCAnt) = real(mFrameTxUp(:,iCAnt) .* exp(1j*2*pi*fCarrFreq*t));
    end

    % ---- F. Output scaling ----
    maxVal = max(abs(mFrameTxCar), [], 'all');
    if maxVal > 0
        mFrameTxCar = 1.5/maxVal * mFrameTxCar;
    end

    % ---- G. Estimate total playback time ----
    wavTime_sec = size(mFrameTxCar,1) / fDACFreq;
    wavTime_ms  = wavTime_sec * 1000 * 1.3;
    wavTime     = ceil(wavTime_ms)/1000;
    meta.wavTime = wavTime;
end