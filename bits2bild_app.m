function bild = bits2bild_app(empfBits, kanal)
%BITS2BILD_APP  Reconstruct image matrix from received bitstream.
%
% Inputs:
%   empfBits : column/row vector of bits (header + image data)
%   kanal    : struct, expected fields
%              .len_cInfoBits  – control info length in bits
%              .bild_size      – [unused, height, width, channels] (optional)
%              .code           – 'Keine' | 'Faltungscodes 2/3' | 'Faltungscodes 1/2'
%
% Output:
%   bild     : uint8 image [H x W x C]

    % --- 0) Basic checks / sanitization ---
    empfBits = empfBits(:);   % force column vector
    nTotal   = numel(empfBits);

    if ~isfield(kanal, 'len_cInfoBits') || isempty(kanal.len_cInfoBits)
        len_cInfoBits = 0;
    else
        len_cInfoBits = kanal.len_cInfoBits;
    end
    len_cInfoBits = max(0, min(len_cInfoBits, nTotal));

    % --- 1) Split control info bits and payload bits ---
    vInfoBits = empfBits(1:len_cInfoBits);
    vEmpfBits = empfBits(len_cInfoBits+1:end);

    % --- 2) Decode control info: image size etc. ---
    bildInfo = [];

    if ~isempty(vInfoBits)
        try
            trel  = poly2trellis(3,[6 7]);   % same code as in bild2bits_app
            tblen = 4;
            decoded_vInfoBits = vitdec(vInfoBits, trel, tblen, 'cont', 'hard');
            decoded_vInfoBits = decoded_vInfoBits(tblen+1:end);

            % Expect 4 fields, each 16 bits
            if numel(decoded_vInfoBits) >= 64
                decoded_vInfoBits = decoded_vInfoBits(1:64);
                bildInfoBits = reshape(decoded_vInfoBits, 16, 4).';
                bildInfo     = bi2de(bildInfoBits).';
            end
        catch
            bildInfo = [];
        end
    end

    % Fallback: use kanal.bild_size if header decode failed
    if (isempty(bildInfo) || numel(bildInfo) ~= 4) ...
            && isfield(kanal,'bild_size') && numel(kanal.bild_size) == 4
        bildInfo = kanal.bild_size(:);
        disp('bits2bild_app: InfoData corrected using kanal.bild_size.');
    end

    % If still no valid size, abort with empty image
    if isempty(bildInfo) || numel(bildInfo) ~= 4 ...
            || any(bildInfo(2:4) <= 0)
        warning('bits2bild_app: unable to recover valid image size, returning empty image.');
        bild = uint8([]);
        return;
    end

    % Number of pixels = H * W * C
    height  = double(bildInfo(2));
    width   = double(bildInfo(3));
    chans   = double(bildInfo(4));
    lgBild  = height * width * chans;

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
            decodedTemp      = vitdec(vEmpfBits, trel, tblen, 'cont', 'hard');
            decoded_vEmpfBits = decodedTemp(tblen+1:end);
        end

    else
        error('bits2bild_app: unknown kanal.code = %s', codeStr);
    end

    % --- 4) Extract image data bits: 8 bits per pixel ---
    neededBits = lgBild * 8;

    if numel(decoded_vEmpfBits) < neededBits
        warning('bits2bild_app: data bits are not enough (%d < %d), zero-padding.', ...
                numel(decoded_vEmpfBits), neededBits);
        vBildBits = [decoded_vEmpfBits(:); zeros(neededBits - numel(decoded_vEmpfBits), 1)];
    else
        vBildBits = decoded_vEmpfBits(1:neededBits);
    end

    % Defensive: ensure multiple of 8
    vBildBits = vBildBits(1:floor(length(vBildBits)/8)*8);
    if isempty(vBildBits)
        bild = uint8([]);
        return;
    end

    % --- 5) 8 bits -> uint8 pixels -> reshape to image ---
    mBildBits = reshape(vBildBits, 8, []).';
    mSymbBild = bi2de(mBildBits);   % [Npix x 1]

    bild = uint8(reshape(mSymbBild, ...
                         height, ...
                         width, ...
                         chans));
end