function [mFrameTxCar, meta] = generateTx_Alamouti_app(params, mDataTxFreq, vPreambleTime)
%GENERATETX_ALAMOUTI_APP
% Alamouti 2x? transmitter:
%   1) Frequency-domain Alamouti mapping (2 time slots per block)
%   2) IFFT to time domain
%   3) Frame structure per original legacy design (5 OFDM symbols per block for Nt=2):
%        col1: Tx1 preamble (only on Tx1)
%        col2: slot1 data (both Tx)
%        col3: slot2 data (both Tx)
%        col4: Tx2 preamble (only on Tx2)
%        col5: null
%   4) Add CP to each OFDM symbol
%   5) Serialize, prepend Schmidl & Cox sync block (no CP), add head/tail zeros
%   6) Resample to DAC rate and upconvert to carrier
%
% Inputs:
%   params      : struct with fields iNoBlocks,iNfft,iNg,iNb,iNoTxAnt,fBBFreq,fDACFreq,fCarrFreq
%   mDataTxFreq : [iNfft x iNoBlocks x 2] frequency-domain data symbols
%   vPreambleTime: [iNfft x 1] time-domain preamble OFDM symbol (no CP)
%
% Outputs:
%   mFrameTxCar : [Nsamples x 2] real passband waveform at fDACFreq
%   meta        : struct with debug info

    % --- Unpack parameters ---
    iNoBlocks = params.iNoBlocks;
    iNfft     = params.iNfft;
    iNg       = params.iNg;

    if isfield(params,'iNb') && ~isempty(params.iNb)
        iNb = params.iNb;
    else
        iNb = iNfft + iNg;
    end

    iNoTxAnt  = params.iNoTxAnt;
    fBBFreq   = params.fBBFreq;
    fDACFreq  = params.fDACFreq;
    fCarrFreq = params.fCarrFreq;

    % --- Minimal checks ---
    assert(iNoTxAnt == 2, 'generateTx_Alamouti_app: iNoTxAnt must be 2 for Alamouti.');
    assert(size(mDataTxFreq,1) == iNfft, 'generateTx_Alamouti_app: mDataTxFreq must have size [iNfft x iNoBlocks x 2].');
    assert(size(mDataTxFreq,2) == iNoBlocks, 'generateTx_Alamouti_app: mDataTxFreq second dim must equal iNoBlocks.');
    assert(size(mDataTxFreq,3) == 2, 'generateTx_Alamouti_app: mDataTxFreq third dim must be 2.');
    assert(numel(vPreambleTime) == iNfft, 'generateTx_Alamouti_app: vPreambleTime length must equal iNfft.');

    vPreambleTime = vPreambleTime(:);

    % --- Defaults to match your Rx timing approach ---
    if ~isfield(params,'headZeros'), params.headZeros = 1000; end
    if ~isfield(params,'tailZeros'), params.tailZeros = 1000; end
    if ~isfield(params,'peak'),      params.peak      = 0.9;  end
    if ~isfield(params,'addSync'),   params.addSync   = true; end

    meta = struct();

    %% 1) Frequency-domain Alamouti mapping
    nBlocksAlam = 2 * iNoBlocks;                 % 2 time slots per original block
    MimoDataTxFreq = zeros(iNfft, nBlocksAlam, 2);

    for iCB = 1:iNoBlocks
        % Slot 1
        MimoDataTxFreq(:, 2*iCB-1, :) = mDataTxFreq(:, iCB, :);

        % Slot 2: [-conj(s2), conj(s1)]
        MimoDataTxFreq(:, 2*iCB, 1) = conj(-mDataTxFreq(:, iCB, 2));
        MimoDataTxFreq(:, 2*iCB, 2) = conj( mDataTxFreq(:, iCB, 1));
    end

    % Power normalization over Tx antennas
    MimoDataTxFreq = (1/sqrt(iNoTxAnt)) * MimoDataTxFreq;

    %% 2) IFFT to time domain (no CP)
    MimoDataTxTime = sqrt(iNfft) * ifft(MimoDataTxFreq, iNfft, 1);  % [iNfft x 2*iNoBlocks x 2]

    %% 3) Frame structure (5 columns per original block for Nt=2)
    nColsPerBlock = iNoTxAnt + 3;   % = 5
    nFrameCols    = nColsPerBlock * iNoBlocks;

    mFrameTxTp = zeros(iNfft, nFrameCols, 2);

    for iCB = 1:iNoBlocks
        vInd = (nColsPerBlock*(iCB-1)+1) : (nColsPerBlock*iCB);  % length 5

        % Insert preambles at legacy positions:
        % Tx1 preamble -> vInd(1)
        % Tx2 preamble -> vInd(2)+2 = vInd(4)
        for iCAnt = 1:2
            colP = vInd(iCAnt) + (iCAnt-1)*2;
            mFrameTxTp(:, colP, iCAnt) = vPreambleTime;
        end

        % Insert two Alamouti slots into columns 2 and 3
        mFrameTxTp(:, vInd(end)-3 : vInd(end)-2, :) = ...
            MimoDataTxTime(:, 2*iCB-1 : 2*iCB, :);

        % Column vInd(end) stays zero (null)
    end

    meta.nColsPerBlock = nColsPerBlock;
    meta.nFrameCols    = nFrameCols;

    %% 4) Add cyclic prefix to each OFDM symbol
    % Result: [iNb x nFrameCols x 2]
    mFrameBB = [mFrameTxTp(end-iNg+1:end,:,:); mFrameTxTp];
    assert(size(mFrameBB,1) == iNb, 'generateTx_Alamouti_app: CP insertion produced unexpected iNb.');

    iNewNoBlocks = size(mFrameBB,2);
    meta.iNewNoBlocks = iNewNoBlocks;

    payloadBB = reshape(mFrameBB, iNb * iNewNoBlocks, 2);   % [NsBB x 2]
    meta.payloadLenBB = size(payloadBB,1);

    %% 5) Prepend Schmidl & Cox sync block (no CP) + head/tail zeros
    mFrameTxBB = payloadBB;

    if params.addSync
        syncNoCp = localBuildSchmidlCoxSync(iNfft);         % [iNfft x 1] complex
        syncBlock = repmat(syncNoCp, 1, 2);                 % [iNfft x 2]
        meta.syncLen = size(syncBlock,1);
    else
        syncBlock = zeros(0,2);
        meta.syncLen = 0;
    end

    headZeros = zeros(params.headZeros, 2);
    tailZeros = zeros(1000, 2);

    mFrameTxBB = [headZeros; syncBlock; mFrameTxBB; tailZeros];

    meta.headLen     = size(headZeros,1);
    meta.tailLen     = size(tailZeros,1);
    meta.basebandLen = size(mFrameTxBB,1);

    %% 6) Resample baseband from fBBFreq -> fDACFreq
    if fBBFreq <= 0 || fDACFreq <= 0
        error('generateTx_Alamouti_app: fBBFreq and fDACFreq must be positive.');
    end

    nBB = size(mFrameTxBB,1);
    nUp = max(1, ceil(nBB * (fDACFreq / fBBFreq)));
    mFrameTxUp = zeros(nUp, 2);

    if exist('resample','file') == 2
        for iCAnt = 1:2
            mFrameTxUp(:,iCAnt) = resample(mFrameTxBB(:,iCAnt), fDACFreq, fBBFreq);
        end
    else
        % Fallback: simple interpolation (runnable but not ideal)
        tBB = (0:nBB-1).' / fBBFreq;
        tUp = (0:nUp-1).' / fDACFreq;
        for iCAnt = 1:2
            mFrameTxUp(:,iCAnt) = interp1(tBB, mFrameTxBB(:,iCAnt), tUp, 'linear', 0);
        end
    end

    meta.upLen = size(mFrameTxUp,1);

    %% 7) Upconvert to carrier frequency (real passband)
    Ns = size(mFrameTxUp,1);
    t  = (0:Ns-1).' / fDACFreq;

    mFrameTxCar = zeros(Ns, 2);
    for iCAnt = 1:2
        mFrameTxCar(:,iCAnt) = real(mFrameTxUp(:,iCAnt) .* exp(1j*2*pi*fCarrFreq*t));
    end

    %% 8) Peak normalization (avoid clipping)
    maxVal = max(abs(mFrameTxCar), [], 'all');
    if maxVal > 0 && isfinite(maxVal)
        mFrameTxCar = (params.peak / maxVal) * mFrameTxCar;
    end

    meta.peakAfterNorm = max(abs(mFrameTxCar), [], 'all');
    meta.wavTimeSec    = size(mFrameTxCar,1) / fDACFreq;
end

% -------------------------------------------------------------------------
% Local Schmidl & Cox sync generator (length = iNfft, no CP)
% A classic choice: two identical halves in time domain.
% -------------------------------------------------------------------------
function syncNoCp = localBuildSchmidlCoxSync(iNfft)
    % Create a random QPSK sequence on half length and repeat it
    L = iNfft/2;
    if mod(iNfft,2) ~= 0
        error('localBuildSchmidlCoxSync: iNfft must be even.');
    end

    % QPSK symbols
    bits = randi([0 1], 2*L, 1);
    sym  = (2*bits(1:2:end)-1) + 1j*(2*bits(2:2:end)-1);
    sym  = sym / sqrt(2);

    % Two identical halves in time
    syncNoCp = [sym; sym];

    % Normalize energy
    syncNoCp = syncNoCp / max(rms(syncNoCp), eps);
end
