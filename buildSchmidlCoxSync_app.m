function [syncSig, Lsync] = buildSchmidlCoxSync_app(iNfft)
%BUILDSCHMIDLCOXSYNC_APP
% Build Schmidl & Cox sync symbol (no CP), length iNfft = 2*L.
% Returns a repeated structure [A; A].

    Lsync = iNfft/2;

    % Complex Gaussian sequence A (length Lsync)
    a = (randn(Lsync,1) + 1j*randn(Lsync,1)) / sqrt(2);

    % Repeated pattern [A; A]
    syncSig = [a; a];

    % Normalize RMS power
    syncSig = syncSig / rms(syncSig);
end