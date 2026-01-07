function [mFrameTxCar, meta] = generateTxSequence_app(params, source)
%GENERATETXSEQUENCE_APP
% Main Tx signal generation:
%   1) file -> bits -> coding -> QAM
%   2) (Text only) overwrite leading symbols with BPSK header
%   3) Chu preamble
%   4) Per-MIMO-mode Tx mapping

    rng(0,'twister');

    % ==== 1. Unpack parameters ====
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iNg          = params.iNg;
    iNb          = params.iNb;
    iModOrd      = params.iModOrd;
    iNoTxAnt     = params.iNoTxAnt;

    mimoModeStr  = string(params.mimoMode);
    coding       = string(params.codingMode);
    datenTyp     = lower(string(source.DatenTyp));

    metaBits = struct();

    % ==== 2. Max bits per frame ====
    anzMaxBits = iModOrd * iNoTxAnt * iNoBlocks * iNfft;

    % ==== 3. File -> bits (with channel coding) ====
    switch datenTyp
        case {"bild","image"}
            [vInfoBits, AnzUeb, metaBild, original_data_bits] = ...
                bild2bits_app(source.SendeDatei, coding, anzMaxBits, iModOrd);
            len_cInfoBits = metaBild.len_cInfoBits;

        case "text"
            [vInfoBits, AnzUeb, metaTxt, original_data_bits] = ...
                text2bits_app(string(source.SendeDatei), coding, anzMaxBits, iModOrd);
            len_cInfoBits = metaTxt.len_cInfoBits;

        otherwise
            error('Unknown DatenTyp: %s (expected Text/Bild/Image).', datenTyp);
    end

    metaBits.AnzUeb           = AnzUeb;
    metaBits.len_cInfoBits    = len_cInfoBits;
    metaBits.originalDataBits = original_data_bits;
    metaBits.anzMaxBits       = anzMaxBits;
    metaBits.DatenTyp         = char(datenTyp);

    % Use a single frame worth of bits
    if length(vInfoBits) < anzMaxBits
        vBits = zeros(1, anzMaxBits);
        vBits(1:length(vInfoBits)) = vInfoBits;
    else
        vBits = vInfoBits(1:anzMaxBits);
    end

    % ==== 4. Serial-to-parallel and symbol indexing ====
    % vBits -> mBits: [iModOrd x (iNoBlocks*iNfft*iNoTxAnt)]
    mBits = reshape(vBits, iModOrd, iNoBlocks*iNfft*iNoTxAnt);

    % IMPORTANT:
    %   Image/Bild: left-msb to match bild2bits_app / bits2bild_app
    %   Text: keep your legacy right-msb behavior to avoid breaking working path
    if datenTyp=="bild" || datenTyp=="image"
        msbOrder = 'left-msb';
    else
        msbOrder = 'right-msb';
    end

    mSymb = bi2de(mBits.', msbOrder);  % column vector indices

    % Frequency-domain symbol indices: [iNfft x iNoBlocks x iNoTxAnt]
    mDataTxFreqTp = reshape(mSymb, iNfft, iNoBlocks, iNoTxAnt);

    % ==== 5. QAM modulation ====
    M = 2^iModOrd;
    mDataTxFreq = zeros(size(mDataTxFreqTp));
    for iCAnt = 1:iNoTxAnt
        mDataTxFreq(:,:,iCAnt) = qammod(mDataTxFreqTp(:,:,iCAnt), M, ...
                                        'UnitAveragePower', true);
    end

    % ==== 6. Control information (header) handling ====
    % Text legacy: overwrite leading symbols with BPSK header (old design)
    % Image/Bild: DO NOT overwrite (header is already embedded in vInfoBits as QAM bits)
    if datenTyp=="text"
        if len_cInfoBits > 0
            vHeaderBits   = vInfoBits(1:len_cInfoBits);
            vDataTxFreq_h = qammod(vHeaderBits, 2, 'UnitAveragePower', true);  % BPSK

            len_einUeb     = iNfft * iNoBlocks * iNoTxAnt;           % symbols per frame
            mDataTxFreqTp2 = reshape(mDataTxFreq, 1, len_einUeb);    % flatten

            if len_cInfoBits <= len_einUeb
                mDataTxFreqTp2(1:len_cInfoBits) = vDataTxFreq_h;
            else
                mDataTxFreqTp2 = vDataTxFreq_h(1:len_einUeb);
            end

            mDataTxFreq = reshape(mDataTxFreqTp2, iNfft, iNoBlocks, iNoTxAnt);
        end
    end

    % ==== 7. Chu preamble in time domain ====
    if exist('ChuSeq','file')
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    else
        vPreambleTime = sqrt(iNoTxAnt) * (randn(iNfft,1) + 1j*randn(iNfft,1))/sqrt(2);
    end

    % ==== 8. Call mode-specific Tx function ====
    modeLower = lower(strtrim(mimoModeStr));

    switch modeLower
        case 'eigenmode'
            [mFrameTxCar, metaMode] = generateTx_EigenMode_app(params, mDataTxFreq, vPreambleTime);

        case 'alamouti'
            [mFrameTxCar, metaMode] = generateTx_Alamouti_app(params, mDataTxFreq, vPreambleTime);

        case {'spatial multiplexing','v-blast','vblast'}
            [mFrameTxCar, metaMode] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime);

        otherwise
            [mFrameTxCar, metaMode] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime);
    end

    meta = mergeMeta(metaBits, metaMode);
end

function metaOut = mergeMeta(metaBase, metaMode)
    metaOut = metaBase;
    if isempty(metaMode), return; end
    fns = fieldnames(metaMode);
    for k = 1:numel(fns)
        fn = fns{k};
        metaOut.(fn) = metaMode.(fn);
    end
end
