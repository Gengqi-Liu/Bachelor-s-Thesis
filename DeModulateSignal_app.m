function vBBSig = DeModulateSignal_app(vHFSig, fBBFreq, fDACFreq, fCarrFreq)
%DEMULATESIGNAL_APP
% RF real signal (vHFSig, fs = fDACFreq) → complex BB signal (fs = fBBFreq)

    % Keep input orientation
    rowInput = isrow(vHFSig);
    vHFSig   = vHFSig(:);

    N = length(vHFSig);

    % 1) Down-conversion
    n        = (0:N-1).';                      % time index
    mix      = exp(-1j*2*pi*fCarrFreq*n/fDACFreq);
    vBBSigUp = vHFSig .* mix;

    % 2) Low-pass filter around baseband bandwidth
    Wn        = min((fBBFreq*1.2) / (fDACFreq/2), 0.99);
    filtOrder = 64;
    B         = fir1(filtOrder, Wn, 'low');

    % 3) Filtering
    vBBSigUpFilt = 2 * filter(B, 1, double(vBBSigUp));

    % 4) Resample to baseband sampling rate
    vBBSig = resample(vBBSigUpFilt, fBBFreq, fDACFreq);

    % Restore original shape
    if rowInput
        vBBSig = vBBSig.';
    end
end