function [bits, AnzUeb, meta, original_data_bits] = bild2bits_app(SendeDatei, codingMode, anzMaxBits, iModOrd)
%BILD2BITS_APP  Image matrix -> bitstream, with optional channel coding.
%
% Inputs:
%   SendeDatei   : image matrix (e.g. output of imread, uint8/uint16)
%   codingMode   : 'Keine' | 'Faltungscodes 2/3' | 'Faltungscodes 1/2'
%   anzMaxBits   : max bits per frame
%   iModOrd      : modulation order (log2(M))
%
% Outputs:
%   bits               : coded + padded bitstream (1 x N)
%   AnzUeb             : number of frames
%   meta               : metadata (image size, control bits etc.)
%   original_data_bits : uncoded bitstream (for BER)

codingMode = string(codingMode);

%% 0. Image -> bitstream
bild = SendeDatei;
[a1, a2, a3] = size(bild);

% Flatten and convert to double for de2bi
vsBild = double(reshape(bild, 1, a1 * a2 * a3));

% Assume 8 bits per pixel component
mBBits = de2bi(vsBild, 8, 'left-msb');
vBBits = reshape(mBBits.', 1, (a1 * a2 * a3) * 8);

original_data_bits   = vBBits;
len_original_data    = length(vBBits);
meta.txt_size        = len_original_data;   % keep naming consistent with text path

%% 1. Channel coding
code_vBits          = vBBits;
Tx_vBits_for_coding = vBBits;

switch codingMode
    case "Keine"
        code_vBits = vBBits;

    case "Faltungscodes 2/3"
        len = length(Tx_vBits_for_coding);

        % Padding and tail bits to make length a multiple of 2
        if mod(len, 2)
            % odd length: +1 data bit +4 tail bits = 5
            Tx_vBits_for_coding = [Tx_vBits_for_coding, randi([0 1], 1, 5)];
        else
            % even length: +4 tail bits
            Tx_vBits_for_coding = [Tx_vBits_for_coding, randi([0 1], 1, 4)];
        end

        % Estimated frame count, stored in first 16 bits
        AnzUeb_est = ceil(length(Tx_vBits_for_coding) / anzMaxBits * 1.5);
        Tx_vBits_for_coding(1:16) = de2bi(AnzUeb_est, 16, 'left-msb');

        % Trellis definition (unchanged from original)
        trel = poly2trellis([4 3], [4 5 17; 7 4 2]);
        code_vBits = convenc(Tx_vBits_for_coding, trel);

    case "Faltungscodes 1/2"
        % Estimated frame count, stored in first 16 bits
        AnzUeb_est = ceil(length(Tx_vBits_for_coding) / anzMaxBits * 2);
        Tx_vBits_for_coding(1:16) = de2bi(AnzUeb_est, 16, 'left-msb');

        trel  = poly2trellis(3, [6 7]);  % K=3, R=1/2
        tblen = 32;                      % tail bits
        tx_vBits_encode = [Tx_vBits_for_coding, randi([0 1], 1, tblen)];
        code_vBits      = convenc(tx_vBits_encode, trel);

    otherwise
        warning('Unknown channel coding mode: %s. Using no coding.', codingMode);
        code_vBits = vBBits;
end

%% 2. Control information: image info (Bild_Info)
% Frame count (based on coded bits length + 136*iModOrd overhead)
AnzUeb = ceil((length(code_vBits) + 136 * iModOrd) / anzMaxBits);

% Image info: [AnzUeb, a1, a2, a3], each 16 bits -> total 64 bits
BildInfo  = de2bi([AnzUeb, a1, a2, a3], 16, 'left-msb');
vbildInfo = reshape(BildInfo.', 1, 64);

meta.bild_size = [AnzUeb, a1, a2, a3];

%% 3. Control information coding (R=1/2 convolutional code)
trel  = poly2trellis(3, [6 7]);  % K=3, R=1/2
tblen = 4;                       % tail bits for control info
bild_infoBits      = [vbildInfo, randi([0 1], 1, tblen)];
bild_infoBits_code = convenc(bild_infoBits, trel);

meta.len_cInfoBits = length(bild_infoBits_code);

%% 4. Prepend control information
% Repeat coded Bild_Info iModOrd times in front of data bits
for i = 1:iModOrd
    code_vBits = [bild_infoBits_code, code_vBits]; %#ok<AGROW>
end

meta.lenCodeBits = length(code_vBits);

%% 5. Pad to AnzUeb * anzMaxBits
len_current      = length(code_vBits);
len_total_needed = AnzUeb * anzMaxBits;
padding_bits     = len_total_needed - len_current;

if padding_bits > 0
    zBits = randi([0 1], 1, padding_bits);
    bits  = [code_vBits, zBits];
else
    bits  = code_vBits;
end
end