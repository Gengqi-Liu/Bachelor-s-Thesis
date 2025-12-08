function [y] = ChuSeq(iLength)
% CHUSEQ Generates a Zadoff-Chu sequence with root r=1 (Chu Sequence).
%
% Input:
% iLength - The length of the sequence (N_ZC).
%
% Output:
% y - The complex Zadoff-Chu sequence (1 x iLength vector).

y = zeros(1,iLength);

% Zadoff-Chu Sequence (r=1) formula: x_u(n) = exp(j * pi * n^2 / N_ZC)
for n = 0:iLength-1
  y(1,n+1) = exp(1j*pi*n*n/iLength);
end
end