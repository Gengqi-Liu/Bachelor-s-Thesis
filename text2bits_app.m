function [bits, AnzUeb, meta, original_data_bits] = text2bits_app(SendeDatei, codingMode, anzMaxBits, iModOrd)
    %TEXT2BITS_APP
    % Text to bitstream with optional convolutional coding.
    
    codingMode = string(codingMode);
    
    % 1) Text -> bits
    symb       = double(char(string(SendeDatei)));
    mBits      = de2bi(symb, 8, 'left-msb');
    len_txBits = numel(mBits);
    tx_vBits   = reshape(mBits.', 1, len_txBits);
    
    original_data_bits = tx_vBits;
    meta               = struct();
    meta.txt_size      = len_txBits;
    
    % 2) Control info (payload length)
    tx_infoBits = de2bi(len_txBits, 64, 'left-msb');
    
    trel_ci  = poly2trellis(3, [6 7]);
    tblen_ci = 4;
    
    tx_infoBits_with_tail = [tx_infoBits, randi([0 1], 1, tblen_ci)];
    tx_infoBits_code      = convenc(tx_infoBits_with_tail, trel_ci);
    
    meta.len_cInfoBits = length(tx_infoBits_code);
    
    % 3) Payload coding
    vBits_to_encode = tx_vBits;
    
    switch codingMode
        case "Keine"
            code_vBits = tx_vBits;
    
        case "Faltungscodes 2/3"
            len_data = length(vBits_to_encode);
    
            if mod(len_data, 2)
                vBits_to_encode = [vBits_to_encode, randi([0 1], 1, 5)];
            else
                vBits_to_encode = [vBits_to_encode, randi([0 1], 1, 4)];
            end
    
            trel_23   = poly2trellis([4 3], [4 5 17; 7 4 2]);
            code_vBits = convenc(vBits_to_encode, trel_23);
    
        case "Faltungscodes 1/2"
            trel_12 = poly2trellis(3, [6 7]);
            tblen   = 32;
    
            vBits_to_encode = [vBits_to_encode, randi([0 1], 1, tblen)];
            code_vBits      = convenc(vBits_to_encode, trel_12);
    
        otherwise
            warning('Unknown coding mode: %s. Using uncoded.', codingMode);
            code_vBits = tx_vBits;
    end
    
    % 4) Prepend control info (repeat iModOrd times)
    for i = 1:iModOrd
        code_vBits = [tx_infoBits_code, code_vBits]; %#ok<AGROW>
    end
    meta.lenCodeBits = length(code_vBits);
    
    % 5) Frame count and padding
    if ~(isscalar(anzMaxBits) && isnumeric(anzMaxBits) && anzMaxBits > 0)
        error('anzMaxBits must be a positive scalar.');
    end
    
    AnzUeb           = ceil(meta.lenCodeBits / anzMaxBits);
    len_total_needed = AnzUeb * anzMaxBits;
    padding_bits     = len_total_needed - meta.lenCodeBits;
    
    if padding_bits > 0
        bits = [code_vBits, randi([0 1], 1, padding_bits)];
    else
        bits = code_vBits;
    end
end