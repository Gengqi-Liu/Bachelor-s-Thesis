function mDataRxEq = VBlastRx_app(mDataRx, mCTF, params, iSNR, kanal)
%VBLASTRX_APP
% V-BLAST equalization (ZF/MMSE). Channel must be provided (mCTF).

    % Unpack
    iNoBlocks    = params.iNoBlocks;
    iNfft        = params.iNfft;
    iNoTxAnt     = params.iNoTxAnt;
    iNoSubBlocks = params.iNoSubBlocks;

    % Validate / clamp
    [Nfft2, nBlocksData, ~] = size(mDataRx);
    if Nfft2 ~= iNfft
        error('VBlastRx_app: Dimension mismatch (iNfft vs mDataRx).');
    end
    if nBlocksData < iNoBlocks
        warning('VBlastRx_app: Fewer data blocks than expected (%d < %d). Using available blocks.', ...
                nBlocksData, iNoBlocks);
        iNoBlocks = nBlocksData;
    end

    % FFT: [SC x Block x Rx]
    mDataRxFreq = 1/sqrt(iNfft) * fft(mDataRx, iNfft, 1);

    % Output: [Tx x SC x Block]
    mDataRxEq = zeros(iNoTxAnt, iNfft, iNoBlocks);

    % SNR guard
    if ~isfinite(iSNR) || iSNR <= 0
        iSNR = 1e4;
    end

    [~, AnzSubFrames, ~, ~] = size(mCTF);

    for iCB = 1:iNoBlocks
        % Subframe index for this block
        iSF = ceil(iCB / iNoSubBlocks);
        if iSF > AnzSubFrames
            iSF = AnzSubFrames;
        end

        for iCSC = 1:iNfft
            % y: [Rx x 1]
            y1 = sqrt(iNoTxAnt) * squeeze(mDataRxFreq(iCSC, iCB, :));

            % H: [Rx x Tx]
            H1 = squeeze(mCTF(iCSC, iSF, :, :));

            % First-stage linear detector
            if strcmpi(kanal.entzerrer, 'Zero Forcing')
                G1 = pinv(H1);
            elseif strcmpi(kanal.entzerrer, 'MMSE')
                A = H1' * H1 + eye(iNoTxAnt) * (1./iSNR);
                if any(~isfinite(A(:))) || rcond(A) < 1e-10
                    A = H1' * H1 + eye(iNoTxAnt) * 1e-3;
                end
                G1 = A \ (H1');
            else
                error('VBlastRx_app: Unknown equalizer: %s', kanal.entzerrer);
            end

            x1 = G1 * y1;  % [Tx x 1]

            % Ordering by row-norm of G1
            vNormSq = sum(abs(G1.').^2, 2).';   % [1 x Tx]
            [~, iFirst] = min(vNormSq);
            [~, iNext]  = max(vNormSq);

            % SIC
            y2 = y1 - H1(:, iFirst) * x1(iFirst);

            H2 = H1;
            H2(:, iFirst) = 0;

            if strcmpi(kanal.entzerrer, 'Zero Forcing')
                G2 = pinv(H2);
            else
                A2 = H2' * H2 + eye(iNoTxAnt) * (1./iSNR);
                if any(~isfinite(A2(:))) || rcond(A2) < 1e-10
                    A2 = H2' * H2 + eye(iNoTxAnt) * 1e-3;
                end
                G2 = A2 \ (H2');
            end

            x2 = G2 * y2;

            % Store
            if iFirst < iNext
                mDataRxEq(iFirst, iCSC, iCB) = x1(iFirst);
                mDataRxEq(iNext,  iCSC, iCB) = x2(iNext);
            else
                mDataRxEq(iNext,  iCSC, iCB) = x2(iNext);
                mDataRxEq(iFirst, iCSC, iCB) = x1(iFirst);
            end
        end
    end
end