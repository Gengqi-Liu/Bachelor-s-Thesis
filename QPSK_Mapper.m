function mDataTxFreq = QPSK_Mapper(vInfoBits, Pglobal)
    %QPSK_MAPPER
    % QPSK (4-QAM) bit-to-symbol mapper.
    %
    % Input:
    %   vInfoBits : bitstream (1 x N) or (N x 1)
    %   Pglobal   : struct with fields .iNfft and .iNoTxAnt
    %
    % Output:
    %   mDataTxFreq : QPSK symbols formatted as [iNfft x iNoBlocks x iNoTxAnt]
    
    % 1) Parameters
    iModOrd  = 2;
    iNfft    = Pglobal.iNfft;
    iNoTxAnt = Pglobal.iNoTxAnt;
    
    vInfoBits = vInfoBits(:).';  % row vector
    
    if mod(numel(vInfoBits), iModOrd) ~= 0
        error('QPSK requires the number of bits to be a multiple of 2.');
    end
    
    % 2) Gray-coded constellation (unit average power)
    norm_factor = 1 / sqrt(2);
    constellation = [ ...
        ( 1+1j) * norm_factor, ... % 00
        ( 1-1j) * norm_factor, ... % 01
        (-1+1j) * norm_factor, ... % 10
        (-1-1j) * norm_factor  ... % 11
    ];
    
    % Bits -> symbol indices
    vBitsReshaped  = reshape(vInfoBits, iModOrd, []).';              % [Nsym x 2]
    vDecimalIdx    = bi2de(vBitsReshaped, 'left-msb') + 1;          % 1..4
    vTxSymbols     = constellation(vDecimalIdx);                    % [1 x Nsym]
    iTotalSymbols  = numel(vTxSymbols);
    
    % 3) Pack into OFDM blocks
    iSymbolsPerBlock = iNfft;
    iNoBlocks        = ceil(iTotalSymbols / iSymbolsPerBlock);
    
    n_pad = iNoBlocks * iSymbolsPerBlock - iTotalSymbols;
    vTxSymbolsPadded = [vTxSymbols, zeros(1, n_pad)];
    
    mSymbols2D = reshape(vTxSymbolsPadded, iSymbolsPerBlock, iNoBlocks);
    
    % 4) Replicate across Tx antennas
    mDataTxFreq = repmat(mSymbols2D, 1, 1, iNoTxAnt);
    
    % Note: do not modify Pglobal here (avoid side effects)
end