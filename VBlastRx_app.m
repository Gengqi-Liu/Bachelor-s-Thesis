function mDataRxEq = VBlastRx_app(mDataRx, mCTF, params, iSNR, kanal)
%VBlastRx_app
% V-BLAST SIC detector per subcarrier.
% Input:
%   mDataRx : [iNfft x iNoBlocks x iNoRxAnt]  (time-domain OFDM blocks, CP removed)
%   mCTF    : [iNfft x AnzSubFrames x iNoRxAnt x iNoTxAnt]
%   params  : struct with fields:
%             iNoBlocks, iNfft, iNoTxAnt, iNoRxAnt, iNoSubBlocks, iModOrd
%   iSNR    : linear SNR estimate (used for MMSE)
%   kanal   : struct with field kanal.entzerrer = 'Zero Forcing'|'MMSE'
%
% Output:
%   mDataRxEq : [iNoTxAnt x iNfft x iNoBlocks] equalized symbols

    % Unpack
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iNoTxAnt     = params.iNoTxAnt;
    iNoRxAnt     = params.iNoRxAnt;
    iNoSubBlocks = params.iNoSubBlocks;

    if isfield(params,'iModOrd') && ~isempty(params.iModOrd)
        iModOrd = params.iModOrd;
    else
        error('VBlastRx_app: params.iModOrd is required for symbol slicing.');
    end
    M = 2^iModOrd;

    % Validate input sizes
    [Nfft2, nBlocksData, Nr2] = size(mDataRx);
    if Nfft2 ~= iNfft
        error('VBlastRx_app: Dimension mismatch (iNfft vs mDataRx).');
    end
    if Nr2 ~= iNoRxAnt
        error('VBlastRx_app: Dimension mismatch (iNoRxAnt vs mDataRx).');
    end
    if nBlocksData < iNoBlocks
        warning('VBlastRx_app: Fewer data blocks than expected (%d < %d). Using available blocks.', ...
                nBlocksData, iNoBlocks);
        iNoBlocks = nBlocksData;
    end

    % SNR guard
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    % Check channel dimensions
    [~, AnzSubFrames, NrH, NtH] = size(mCTF);
    if NrH ~= iNoRxAnt
        error('VBlastRx_app: Dimension mismatch (mCTF Rx dim %d vs iNoRxAnt %d).', NrH, iNoRxAnt);
    end
    if NtH ~= iNoTxAnt
        error('VBlastRx_app: Dimension mismatch (mCTF Tx dim %d vs iNoTxAnt %d).', NtH, iNoTxAnt);
    end

    % FFT: [SC x Block x Rx]
    mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);

    % Output: [Tx x SC x Block]
    mDataRxEq = zeros(iNoTxAnt, iNfft, iNoBlocks);

    for iCB = 1:iNoBlocks
        % Subframe index for this block
        iSF = ceil(iCB / iNoSubBlocks);
        if iSF > AnzSubFrames
            iSF = AnzSubFrames;
        end

        for iCSC = 1:iNfft
            % y: [Rx x 1]
            y = sqrt(iNoTxAnt) * squeeze(mDataRxFreq(iCSC, iCB, :));

            % H: [Rx x Tx]
            H = squeeze(mCTF(iCSC, iSF, :, :));

            % Full SIC
            Hrem = H;
            yrem = y;
            idxRemain = 1:iNoTxAnt;

            xHatFull = zeros(iNoTxAnt,1);

            for it = 1:iNoTxAnt
                % Linear filter on remaining system
                if strcmpi(kanal.entzerrer, 'Zero Forcing')
                    G = pinv(Hrem);
                elseif strcmpi(kanal.entzerrer, 'MMSE')
                    NtRem = size(Hrem,2);
                    A = Hrem' * Hrem + eye(NtRem) * (1./iSNR);
                    if any(~isfinite(A(:))) || rcond(A) < 1e-10
                        A = Hrem' * Hrem + eye(NtRem) * 1e-3;
                    end
                    G = A \ (Hrem');
                else
                    error('VBlastRx_app: Unknown equalizer: %s', kanal.entzerrer);
                end

                xSoft = G * yrem;  % [NtRem x 1]

                % Ordering: pick stream with minimum post-filter noise enhancement
                % Row-norm of G corresponds to noise enhancement per stream
                vNormSq = sum(abs(G).^2, 2);      % [NtRem x 1]
                [~, posPick] = min(vNormSq);        % position in remaining set
                iPick = idxRemain(posPick);         % original stream index

                % Hard slicing to nearest constellation point
                symIdx = qamdemod(xSoft(posPick), M, 'UnitAveragePower', true);
                xHard  = qammod(symIdx, M, 'UnitAveragePower', true);

                xHatFull(iPick) = xHard;

                % Cancel and remove detected stream
                yrem = yrem - Hrem(:,posPick) * xHard;
                Hrem(:,posPick) = [];
                idxRemain(posPick) = [];

                if isempty(idxRemain)
                    break;
                end
            end

            mDataRxEq(:, iCSC, iCB) = xHatFull;
        end
    end
end
