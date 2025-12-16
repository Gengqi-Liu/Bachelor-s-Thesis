function [frameStart, fCarOffset, MetricData] = EstFrameStart_app(mDataRx, iNfft, iNg)
% EstFrameStart_app
%
% Schmidl & Cox timing sync + coarse CFO estimate.
%
% Input:
%   mDataRx : 1xN complex baseband (single channel)
%   iNfft   : sync symbol length (= 2L, must be even)
%   iNg     : unused (kept for compatibility)
%
% Output:
%   frameStart : estimated start index of the sync symbol (1-based)
%   fCarOffset : coarse CFO (normalized frequency, cycles/sample)
%   MetricData : struct with P, R, Metric, idxMax

    %#ok<INUSD>
    r = mDataRx(:).';
    Nsig = numel(r);

    if mod(iNfft, 2) ~= 0
        error('EstFrameStart_app: iNfft must be even (iNfft = 2*L).');
    end
    L = iNfft/2;

    maxD = Nsig - 2*L;
    if maxD <= 0
        error('EstFrameStart_app: input too short for one full sync symbol.');
    end

    P = zeros(1, maxD);
    R = zeros(1, maxD);

    for d = 1:maxD
        seg1 = r(d     : d+L-1);
        seg2 = r(d+L   : d+2*L-1);

        P(d) = sum(conj(seg1) .* seg2);
        R(d) = sum(abs(seg1).^2 + abs(seg2).^2);
    end

    Metric = (abs(P).^2) ./ (R.^2 + eps);

    [~, idxMax] = max(Metric);
    frameStart  = idxMax;

    fCarOffset = angle(P(idxMax)) / (2*pi*L);

    MetricData = struct();
    MetricData.P      = P;
    MetricData.R      = R;
    MetricData.Metric = Metric;
    MetricData.idxMax = idxMax;
end