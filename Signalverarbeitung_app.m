function data = Signalverarbeitung_app(procParam)
% Signalverarbeitung_app
%
% End-to-end RX processing pipeline (baseband conversion, sync, channel
% estimation, equalization, bit recovery, and payload reconstruction).
%
% Input:
%   procParam struct with fields:
%     rxSignal          [Nsamp x Nr] recorded waveform (audioPlayerRecorder)
%     fs                sample rate (= fDACFreq)
%     iNfft, iNg, iNb
%     iNoBlocks, iNewNoBlocks, iNoSubBlocks
%     iNoTxAnt, iNoRxAnt
%     iModOrd
%     fBBFreq, fCarrFreq
%     mimoMode
%     codingMode
%     channelEstimator
%     equalizerMode
%     len_cInfoBits     (optional) control/header bit length (default: 0)
%     DatenTyp          'Bild'/'Image'/'Text'
%     SendeDatei        original payload (image matrix or text)
%
% Output:
%   data struct containing:
%     text / bild
%     numBitError, bitErrorRate
%     plus intermediate PHY results returned by AnalyzeRxSig_app

    %% 1) Unpack basic inputs
    rxSignal = procParam.rxSignal;   % [Nsamp x Nr]
    fs       = procParam.fs;

    %% 2) Build parameter structs for PHY processing
    params = struct( ...
        'iNoBlocks',    procParam.iNoBlocks, ...
        'iNfft',        procParam.iNfft, ...
        'iNg',          procParam.iNg, ...
        'iNb',          procParam.iNb, ...
        'iModOrd',      procParam.iModOrd, ...
        'iNoTxAnt',     procParam.iNoTxAnt, ...
        'iNoRxAnt',     procParam.iNoRxAnt, ...
        'iNewNoBlocks', procParam.iNewNoBlocks, ...
        'iNoSubBlocks', procParam.iNoSubBlocks, ...
        'fBBFreq',      procParam.fBBFreq, ...
        'fDACFreq',     fs, ...
        'fCarrFreq',    procParam.fCarrFreq, ...
        'mimoMode',     procParam.mimoMode ...
    );

    if isfield(procParam,'len_cInfoBits')
        len_cInfoBits = procParam.len_cInfoBits;
    else
        len_cInfoBits = 0;
    end

    kanal = struct( ...
        'len_cInfoBits', len_cInfoBits, ...
        'code',          procParam.codingMode, ...
        'schaetzer',     procParam.channelEstimator, ...
        'entzerrer',     procParam.equalizerMode ...
    );

    %% 3) PHY processing: RF -> BB -> sync -> channel -> equalization -> bits
    % AnalyzeRxSig_app expects rxFrame as [Nr x Nsamp]
    [mEmpfDataBits, dataChan] = AnalyzeRxSig_app(rxSignal.', params, kanal);

    data = dataChan;

    %% 4) Bits -> payload reconstruction (text/image) + BER
    dt = lower(string(procParam.DatenTyp));

    switch dt
        case {"bild","image"}
            bild = bits2bild_app(mEmpfDataBits, kanal);
            data.bild = bild;
            data.text = '';

            if isfield(procParam,'SendeDatei') && ~isempty(procParam.SendeDatei)
                [data.numBitError, data.bitErrorRate] = ...
                    bitFehlerRaten_app(bild, procParam.SendeDatei, procParam.DatenTyp);
            else
                data.numBitError  = NaN;
                data.bitErrorRate = NaN;
            end

        case "text"
            [textDecoded, ~] = bits2text_app(mEmpfDataBits, kanal);
            data.text = textDecoded;
            data.bild = [];

            if isfield(procParam,'SendeDatei') && ~isempty(procParam.SendeDatei)
                [data.numBitError, data.bitErrorRate] = ...
                    bitFehlerRaten_app(textDecoded, procParam.SendeDatei, procParam.DatenTyp);
            else
                data.numBitError  = NaN;
                data.bitErrorRate = NaN;
            end

        otherwise
            data.text = '(Unknown DatenTyp: cannot reconstruct payload)';
            data.bild = [];
            data.numBitError  = NaN;
            data.bitErrorRate = NaN;
    end
end