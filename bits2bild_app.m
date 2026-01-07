function bild = bits2bild_app(empfBits, kanal)
%BITS2BILD_APP  Reconstruct image matrix from received bitstream.
%
% Inputs:
%   empfBits : column/row vector of bits (header + image data)
%   kanal    : struct, expected fields
%              .len_cInfoBits  – control info length in bits (for ONE copy)
%              .code           – 'Keine' | 'Faltungscodes 2/3' | 'Faltungscodes 1/2'
%              .iModOrd        – optional, bits per QAM symbol (1..6)
%
% Output:
%   bild     : uint8 image [H x W x C]

    % --- 0) Basic checks / sanitization ---
    empfBits = empfBits(:);   % force column
    nTotal   = numel(empfBits);

    if nTotal == 0
        bild = uint8([]);
        return;
    end

    if ~isfield(kanal, 'len_cInfoBits') || isempty(kanal.len_cInfoBits)
        len_cInfoBits = 0;
    else
        len_cInfoBits = kanal.len_cInfoBits;
    end
    len_cInfoBits = max(0, min(len_cInfoBits, nTotal));

    % Optional: iModOrd (header repeated iModOrd times on TX)
    if isfield(kanal,'iModOrd') && ~isempty(kanal.iModOrd) && isfinite(kanal.iModOrd)
        iModOrd = round(double(kanal.iModOrd));
        iModOrd = min(max(iModOrd,1),6);
    else
        iModOrd = 1; % safe default if not provided
    end

    % --- Helper: decode one header chunk (len_cInfoBits bits) -> [AnzUeb,H,W,C] ---
    function info = decodeHeaderChunk(vInfoBitsChunk)
        info = [];
        if isempty(vInfoBitsChunk)
            return;
        end
        try
            trel  = poly2trellis(3,[6 7]);   % same as TX control code
            tblen = 4;                       % same tail length as TX for control

            % Use truncation for a finite-length block
            dec = vitdec(vInfoBitsChunk(:).', trel, tblen, 'trunc', 'hard');
            dec = dec(:);

            % Drop tail bits (TX appended 4 tail bits before coding)
            if numel(dec) > tblen
                dec = dec(tblen+1:end);
            else
                return;
            end

            if numel(dec) < 64
                return;
            end

            dec = dec(1:64);
            m16 = reshape(dec, 16, 4).';                     % 4 words, 16 bits each
            info = bi2de(m16, 'left-msb').';                 % IMPORTANT: left-msb matches TX
        catch
            info = [];
        end
    end

    % --- 1) Try to decode control info (handle possible repetitions) ---
    bildInfo = [];
    repUsed  = 1;

    if len_cInfoBits >= 1 && nTotal >= len_cInfoBits
        maxRepTry = min(iModOrd, floor(nTotal / len_cInfoBits));
        maxRepTry = min(maxRepTry, 6);  % clamp search

        for r = 1:maxRepTry
            chunk = empfBits((r-1)*len_cInfoBits + (1:len_cInfoBits));
            cand  = decodeHeaderChunk(chunk);

            % Validate candidate header
            if ~isempty(cand) && numel(cand) == 4
                H = double(cand(2));
                W = double(cand(3));
                C = double(cand(4));

                % Hard sanity limits (adjust if you transmit larger images)
                if H >= 1 && H <= 2048 && W >= 1 && W <= 2048 && (C == 1 || C == 3)
                    bildInfo = cand(:);
                    repUsed  = r;
                    break;
                end
            end
        end
    end

    % If header still not recovered, try kanal.bild_size as fallback
    if (isempty(bildInfo) || numel(bildInfo) ~= 4) && isfield(kanal,'bild_size') && numel(kanal.bild_size) == 4
        bildInfo = kanal.bild_size(:);
        repUsed  = iModOrd; % best effort: assume TX repeated it iModOrd times
    end

    if isempty(bildInfo) || numel(bildInfo) ~= 4
        warning('bits2bild_app: Unable to recover a valid image header.');
        bild = uint8([]);
        return;
    end

    height = double(bildInfo(2));
    width  = double(bildInfo(3));
    chans  = double(bildInfo(4));

    if ~(height >= 1 && height <= 2048 && width >= 1 && width <= 2048 && (chans == 1 || chans == 3))
        warning('bits2bild_app: Invalid image size decoded: %dx%dx%d.', height, width, chans);
        bild = uint8([]);
        return;
    end

    % --- 2) Split payload bits (skip ALL header copies used) ---
    hdrTotalBits = repUsed * len_cInfoBits;
    hdrTotalBits = min(hdrTotalBits, nTotal);

    vEmpfBits = empfBits(hdrTotalBits+1:end);

    % --- 3) Channel decoding for image payload bits ---
    if ~isfield(kanal,'code') || isempty(kanal.code)
        codeStr = 'Keine';
    else
        codeStr = char(kanal.code);
    end

    if strcmpi('Keine', codeStr)
        decoded_vEmpfBits = vEmpfBits(:);

    elseif strcmpi('Faltungscodes 2/3', codeStr)
        vEmpfBits = vEmpfBits(:);
        len = floor(length(vEmpfBits)/1.5);
        if len <= 0
            decoded_vEmpfBits = zeros(0,1);
        else
            if mod(len, 2)
                vEmpfBits = vEmpfBits(1:floor((len-1)*1.5));
            else
                vEmpfBits = vEmpfBits(1:floor(len*1.5));
            end
            trel  = poly2trellis([4 3],[4 5 17;7 4 2]);
            tblen = 32;
            decoded_vEmpfBits = vitdec(vEmpfBits, trel, tblen, 'trunc', 'hard');
            decoded_vEmpfBits = decoded_vEmpfBits(:);
        end

    elseif strcmpi('Faltungscodes 1/2', codeStr)
        vEmpfBits = vEmpfBits(:);
        len = length(vEmpfBits);
        if len < 2
            decoded_vEmpfBits = zeros(0,1);
        else
            if mod(len, 2)
                vEmpfBits = vEmpfBits(1:(len-1));
            end
            trel  = poly2trellis(3,[6 7]);
            tblen = 16;
            decodedTemp       = vitdec(vEmpfBits, trel, tblen, 'trunc', 'hard');
            decoded_vEmpfBits = decodedTemp(:);
            if numel(decoded_vEmpfBits) > tblen
                decoded_vEmpfBits = decoded_vEmpfBits(tblen+1:end);
            end
        end

    else
        error('bits2bild_app: Unknown kanal.code = %s', codeStr);
    end

    % --- 4) Extract image data bits: 8 bits per pixel ---
    lgBild     = height * width * chans;
    neededBits = lgBild * 8;

    % Absolute safety guard to prevent huge allocations
    maxBitsAllowed = 50e6; % 50 Mbits
    if neededBits > maxBitsAllowed
        error('bits2bild_app: neededBits=%g is too large. Header is likely corrupted.', neededBits);
    end

    if numel(decoded_vEmpfBits) < neededBits
        % Pad only to a reasonable size
        missing = neededBits - numel(decoded_vEmpfBits);
        if missing > maxBitsAllowed
            error('bits2bild_app: Padding length=%g is too large. Header is likely corrupted.', missing);
        end
        vBildBits = [decoded_vEmpfBits(:); zeros(missing, 1)];
    else
        vBildBits = decoded_vEmpfBits(1:neededBits);
    end

    % Ensure multiple of 8
    vBildBits = vBildBits(1:floor(numel(vBildBits)/8)*8);
    if isempty(vBildBits)
        bild = uint8([]);
        return;
    end

    % --- 5) 8 bits -> uint8 pixels -> reshape ---
    mBildBits = reshape(vBildBits, 8, []).';
    mSymbBild = bi2de(mBildBits, 'left-msb');  % IMPORTANT: left-msb matches TX de2bi
    bild      = uint8(reshape(mSymbBild, height, width, chans));
end
