function [mFrameTxCar, meta] = generateTxSequence_app(params, source)
%GENERATETXSEQUENCE_APP
% Main Tx signal generation:
%   1) Common pipeline: file -> bits -> coding -> QAM -> control header -> Chu preamble
%   2) Per-MIMO-mode Tx mapping (EigenMode / SM / V-BLAST / Alamouti)

    % Fixed random seed (for reproducibility)
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
    datenTyp     = string(source.DatenTyp);

    metaBits = struct();   % bit-level meta information

    % ==== 2. Max bits per frame ====
    anzMaxBits = iModOrd * iNoTxAnt * iNoBlocks * iNfft;

    % ==== 3. File -> bits (with channel coding) ====
    switch lower(datenTyp)
        case "bild"
            if exist('bild2bits_app','file')
                [vInfoBits, AnzUeb, metaBild, original_data_bits] = ...
                    bild2bits_app(source.SendeDatei, coding, anzMaxBits, iModOrd);
            else
                error('bild2bits_app not found. Please add it to the path.');
            end
            len_cInfoBits = metaBild.len_cInfoBits;

        case "image"
            if exist('bild2bits_app','file')
                [vInfoBits, AnzUeb, metaBild, original_data_bits] = ...
                    bild2bits_app(source.SendeDatei, coding, anzMaxBits, iModOrd);
            else
                error('bild2bits_app not found. Please add it to the path.');
            end
            len_cInfoBits = metaBild.len_cInfoBits;

        case "text"
            if exist('text2bits_app','file')
                [vInfoBits, AnzUeb, metaTxt, original_data_bits] = ...
                    text2bits_app(string(source.SendeDatei), coding, anzMaxBits, iModOrd);
            else
                error('text2bits_app not found. Please add it to the path.');
            end
            len_cInfoBits = metaTxt.len_cInfoBits;

        otherwise
            error('Unknown DatenTyp: %s (expected Text/Bild/Image).', datenTyp);
    end

    % Common meta
    metaBits.AnzUeb           = AnzUeb;
    metaBits.len_cInfoBits    = len_cInfoBits;
    metaBits.originalDataBits = original_data_bits;
    metaBits.anzMaxBits       = anzMaxBits;

    % Use a single frame worth of bits
    if length(vInfoBits) < anzMaxBits
        vBits = zeros(1, anzMaxBits);
        vBits(1:length(vInfoBits)) = vInfoBits;
    else
        vBits = vInfoBits(1:anzMaxBits);
    end

    % ==== 4. Serial-to-parallel and symbol indexing ====
    % vBits -> mBits: [iModOrd x (iNoBlocks*iNfft*iNoTxAnt)]
    mBits         = reshape(vBits, iModOrd, iNoBlocks*iNfft*iNoTxAnt);
    % Each column is one symbol (right-msb)
    mSymb         = bi2de(mBits.', 'right-msb');   % column vector
    % Frequency-domain symbol indices: [iNfft x iNoBlocks x iNoTxAnt]
    mDataTxFreqTp = reshape(mSymb, iNfft, iNoBlocks, iNoTxAnt);

    % ==== 5. QAM modulation ====
    M = 2^iModOrd;
    mDataTxFreq = zeros(size(mDataTxFreqTp));
    for iCAnt = 1:iNoTxAnt
        mDataTxFreq(:,:,iCAnt) = qammod(mDataTxFreqTp(:,:,iCAnt), M, ...
                                        'UnitAveragePower', true);
    end

    % ==== 6. Control information (header) in BPSK on leading symbols ====
    vHeaderBits   = vInfoBits(1:len_cInfoBits);
    vDataTxFreq_h = qammod(vHeaderBits, 2, 'UnitAveragePower', true);  % BPSK

    len_einUeb     = iNfft * iNoBlocks * iNoTxAnt;           % symbols per frame
    mDataTxFreqTp2 = reshape(mDataTxFreq, 1, len_einUeb);    % flatten

    if len_cInfoBits <= len_einUeb
        % Header occupies first len_cInfoBits positions
        mDataTxFreqTp2(1:len_cInfoBits) = vDataTxFreq_h;
        anzMaxUeb_h = 1;
    else
        % Header longer than single frame, truncate
        mDataTxFreqTp2 = vDataTxFreq_h(1:len_einUeb);
        anzMaxUeb_h    = ceil(len_cInfoBits / len_einUeb);
    end

    mDataTxFreq          = reshape(mDataTxFreqTp2, iNfft, iNoBlocks, iNoTxAnt);
    metaBits.anzMaxUeb_h = anzMaxUeb_h;

    % ==== 7. Chu preamble in time domain ====
    if exist('ChuSeq','file')
        vPreambleTime = sqrt(iNoTxAnt) * ChuSeq(iNfft).';
    else
        warning('ChuSeq not found. Using random preamble.');
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
            warning('Unknown MIMO mode "%s". Falling back to Spatial Multiplexing.', mimoModeStr);
            [mFrameTxCar, metaMode] = generateTx_SM_app(params, mDataTxFreq, vPreambleTime);
    end

    % ==== 9. Merge meta and return ====
    meta = mergeMeta(metaBits, metaMode);

end

% Helper: merge meta structures
function metaOut = mergeMeta(metaBase, metaMode)
    metaOut = metaBase;
    if isempty(metaMode)
        return;
    end
    fns = fieldnames(metaMode);
    for k = 1:numel(fns)
        fn = fns{k};
        metaOut.(fn) = metaMode.(fn);
    end
end