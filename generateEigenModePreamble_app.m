function [txPreCar, metaPre] = generateEigenModePreamble_app(params)
% generateEigenModePreamble_app
%
% Build an EigenMode channel-sounding preamble and upconvert to carrier.
%
% Input (params):
%   iNfft, iNg
%   iNb        (optional; default iNfft+iNg)
%   iNoTxAnt
%   fBBFreq, fDACFreq, fCarrFreq
%   Npre       (optional; number of repeats, default 1)
%
% Output:
%   txPreCar : [Nsamp x Nt] real-valued passband signal
%   metaPre  : struct with fields:
%       N, iNewNoBlocksPre, headLen, tailLen, bbLen, wavTime

    % ---- Unpack params ----
    iNfft = params.iNfft;
    iNg   = params.iNg;

    if isfield(params,'iNb') && ~isempty(params.iNb)
        iNb = params.iNb;
    else
        iNb = iNfft + iNg;
    end

    iNoTxAnt  = params.iNoTxAnt;
    fBBFreq   = params.fBBFreq;
    fDACFreq  = params.fDACFreq;
    fCarrFreq = params.fCarrFreq;

    if isfield(params,'Npre') && ~isempty(params.Npre)
        N = params.Npre;
    else
        N = 1;
    end

    metaPre = struct();
    metaPre.N = N;

    % ---- Chu preamble in time domain ----
    if exist('ChuSeq','file')
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft);
        vPreambleTime = vPreambleTime(:).';
    else
        warning('generateEigenModePreamble_app: ChuSeq not found, using random sequence.');
        vPreambleTime = sqrt(iNoTxAnt) * (randn(1,iNfft) + 1j*randn(1,iNfft)) / sqrt(2);
    end

    % ---- Per-antenna time-division preamble (legacy style) ----
    preLen = numel(vPreambleTime);
    PreambleFrame = zeros(iNoTxAnt, preLen * iNoTxAnt);

    k = 0;
    for iC = 1:iNoTxAnt
        PreambleFrame(iC, k+1 : k+preLen) = vPreambleTime;
        k = k + preLen;
    end

    % ---- Repeat N times for averaging ----
    PreambleFrame_lang = zeros(iNoTxAnt, preLen * iNoTxAnt * N);
    for iRep = 1:N
        idx = (iRep-1)*preLen*iNoTxAnt + (1:preLen*iNoTxAnt);
        PreambleFrame_lang(:, idx) = PreambleFrame;
    end

    % ---- Add CP per OFDM block (baseband) ----
    iNewNoBlocksPre = N * iNoTxAnt;
    metaPre.iNewNoBlocksPre = iNewNoBlocksPre;

    PreambleFrameCp = zeros(iNb, iNewNoBlocksPre, iNoTxAnt);

    for iC = 1:iNoTxAnt
        tmp = reshape(PreambleFrame_lang(iC,:), iNfft, iNewNoBlocksPre);
        PreambleFrameCp(:,:,iC) = [tmp(end-iNg+1:end,:); tmp];
    end

    mPreambleBB = reshape(PreambleFrameCp, iNb*iNewNoBlocksPre, iNoTxAnt);

    % ---- Head/tail silence in baseband ----
    headZeros = zeros(1000, iNoTxAnt);
    tailZeros = zeros(3000, iNoTxAnt);
    mPreambleBB = [headZeros; mPreambleBB; tailZeros];

    metaPre.headLen = size(headZeros,1);
    metaPre.tailLen = size(tailZeros,1);
    metaPre.bbLen   = size(mPreambleBB,1);

    % ---- Resample to DAC rate ----
    Ns_up  = ceil(size(mPreambleBB,1) * fDACFreq / fBBFreq);
    txPreUp = zeros(Ns_up, iNoTxAnt);

    for iC = 1:iNoTxAnt
        txPreUp(:,iC) = resample(mPreambleBB(:,iC), fDACFreq, fBBFreq);
    end

    % ---- Upconvert to carrier ----
    Ns = size(txPreUp,1);
    t  = (0:Ns-1).' / fDACFreq;

    txPreCar = zeros(Ns, iNoTxAnt);
    for iC = 1:iNoTxAnt
        txPreCar(:,iC) = real(txPreUp(:,iC) .* exp(1j*2*pi*fCarrFreq*t));
    end

    % ---- Peak normalization ----
    maxVal = max(abs(txPreCar), [], 'all');
    if maxVal > 0
        txPreCar = 0.6 / maxVal * txPreCar;
    end

    % ---- Playback time estimate ----
    wavTime_sec = size(txPreCar,1) / fDACFreq;
    wavTime_ms  = wavTime_sec * 1000 * 1.3;
    metaPre.wavTime = ceil(wavTime_ms) / 1000;
end