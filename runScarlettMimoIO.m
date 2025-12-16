function [rx_raw, usedRealDevice] = runScarlettMimoIO( ...
        txSig, fs, Nt, Nr, deviceName, frameLen)
% runScarlettMimoIO
%
% txSig      : [Nsamples x Nt] transmit signal (real or complex)
% fs         : sample rate (Hz), must match audio device
% Nt         : number of Tx channels
% Nr         : number of Rx channels
% deviceName : audio device name (e.g. 'Scarlett 18i20 USB')
% frameLen   : audio I/O frame length (samples)
%
% Outputs:
% rx_raw         : [Nrecv x Nr] received waveform
% usedRealDevice : true  -> real audio device used
%                  false -> fallback to ideal loopback (rx = tx)
%
% Side effects:
%   Saves audio files in current directory:
%       tx_last.wav  - transmitted multichannel audio
%       rx_last.wav  - received multichannel audio

    %% 0) Ensure real-valued signal and apply amplitude protection
    tx = real(txSig);                     % discard imaginary part if present
    maxVal = max(abs(tx), [], 'all');
    if maxVal > 1
        tx = tx ./ maxVal;                % prevent clipping
    end
    txGain = 0.7;                         % fixed linear gain
    tx = txGain * tx;

    Nsamples = size(tx,1);

    %% 1) Try to open real audio device
    rx_raw = [];
    usedRealDevice = false;

    try
        aPR = audioPlayerRecorder( ...
            'Device', deviceName, ...
            'SampleRate', fs, ...
            'PlayerChannelMapping', 1:Nt, ...
            'RecorderChannelMapping', 1:Nr);
        usedRealDevice = true;
        fprintf('Audio device "%s" opened. Starting I/O...\n', deviceName);
    catch
        usedRealDevice = false;
    end

    %% 2) Real-device path: frame-based playback and recording
    if usedRealDevice
        ptr = 1;
        rx_buf = zeros(0, Nr);

        while ptr <= Nsamples
            idxEnd = min(ptr+frameLen-1, Nsamples);
            frame = tx(ptr:idxEnd, :);      % [L x Nt]

            if size(frame,1) < frameLen
                frame(end+1:frameLen,:) = 0; %#ok<AGROW>
            end

            [rframe, underrun, overrun] = aPR(frame);

            if underrun > 0
                warning('Audio underrun = %d', underrun);
            end
            if overrun > 0
                warning('Audio overrun = %d', overrun);
            end

            rx_buf = [rx_buf; rframe]; %#ok<AGROW>
            ptr = ptr + frameLen;
        end

        if exist('aPR','var') && isa(aPR,'audioPlayerRecorder')
            release(aPR);
        end

        rx_raw = rx_buf;
        fprintf('Audio I/O finished: %d samples x %d channels.\n', ...
                size(rx_raw,1), size(rx_raw,2));

        % Save Tx/Rx audio
        try
            audiowrite('tx_last.wav', tx, fs);
            audiowrite('rx_last.wav', rx_raw, fs);
            fprintf('Saved tx_last.wav and rx_last.wav.\n');
        catch MEw
            warning('Failed to save wav files: %s', MEw.message);
        end

    %% 3) Simulation path: ideal noiseless loopback
    else
        ch = min(Nt, Nr);
        rx_raw = zeros(size(tx,1), Nr);
        rx_raw(:,1:ch) = tx(:,1:ch);

        fprintf('Simulation mode: ideal loopback (rx = tx).\n');

        try
            audiowrite('tx_last.wav', tx, fs);
            audiowrite('rx_last.wav', rx_raw, fs);
            fprintf('Saved tx_last.wav and rx_last.wav (simulation mode).\n');
        catch MEw
            warning('Failed to save wav files: %s', MEw.message);
        end
    end
end