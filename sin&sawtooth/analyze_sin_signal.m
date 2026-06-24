% Offline analysis and recovery for the sine signal saved by pureRX.slx.
% This script reads D:\drone\0622\0622\SINsignal.mat and does not modify it.

clear; clc; close all;

matFile = locate_mat_file('SINsignal2.mat');
Fs = 1e6;

rxCube = load_first_signal(matFile);
rx = rxCube(:);
rx = rx(isfinite(real(rx)) & isfinite(imag(rx)));
rx = rx - mean(rx);

fprintf('Sine signal offline analysis\n');
fprintf('  File: %s\n', matFile);
fprintf('  Total samples: %d\n', numel(rx));
fprintf('  Duration: %.3f s\n', numel(rx)/Fs);
fprintf('  RMS amplitude: %.4g\n', rms(abs(rx)));
fprintf('  Peak amplitude: %.4g\n', max(abs(rx)));

[toneFreq, tonePowerDb] = estimate_tone_frequency(rx, Fs);
fprintf('  Estimated received tone offset: %.3f Hz\n', toneFreq);
fprintf('  Tone PSD peak: %.2f dB\n', tonePowerDb);

numPlot = min(numel(rx), 20000);
tPlot = (0:numPlot-1).' / Fs;

figure('Name', 'SINsignal Raw IQ Time Domain');
subplot(3,1,1);
plot(tPlot, real(rx(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('I');
title('Raw received sine signal, real part');

subplot(3,1,2);
plot(tPlot, imag(rx(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('Q');
title('Raw received sine signal, imaginary part');

subplot(3,1,3);
plot(tPlot, abs(rx(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('|r[n]|');
title('Raw received sine signal amplitude');

figure('Name', 'SINsignal Spectrum');
nfftWelch = 8192;
[pxx, f] = pwelch(rx, hamming(nfftWelch), nfftWelch/2, nfftWelch, Fs, 'centered');
plot(f/1e3, 10*log10(pxx + eps));
grid on;
xlabel('Frequency offset (kHz)'); ylabel('PSD (dB/Hz)');
title(sprintf('Received spectrum, peak offset %.3f Hz', toneFreq));

% Coherent tone recovery: mix the estimated tone to DC, smooth its complex
% amplitude, then remodulate to recover a clean sine at the measured offset.
n = (0:numel(rx)-1).';
basebandTone = rx .* exp(-1j*2*pi*toneFreq/Fs*n);

avgLen = 2000;
toneEnvelope = movmean(basebandTone, avgLen);
complexAmp = median(toneEnvelope(round(0.1*numel(toneEnvelope)):end));
recoveredComplex = complexAmp .* exp(1j*2*pi*toneFreq/Fs*n);
recoveredSine = real(recoveredComplex);

% A second output useful for analysis: the slowly varying demodulated
% envelope after the tone is shifted to DC.
demodEnvelope = toneEnvelope;

fprintf('  Recovered sine amplitude estimate: %.4g\n', abs(complexAmp));
fprintf('  Recovered sine phase estimate: %.3f rad\n', angle(complexAmp));

figure('Name', 'Recovered Sine Comparison');
subplot(2,1,1);
plot(tPlot, real(rx(1:numPlot)), 'Color', [0.6 0.6 0.6]);
hold on;
plot(tPlot, recoveredSine(1:numPlot), 'r', 'LineWidth', 1.1);
grid on;
xlabel('Time (s)'); ylabel('Amplitude');
legend('Raw I', 'Recovered sine');
title('Recovered sine over raw I channel');

subplot(2,1,2);
plot(tPlot, abs(demodEnvelope(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('Envelope');
title('Demodulated tone envelope after frequency shift to DC');

figure('Name', 'Recovered Sine Spectrum');
[pxxRec, fRec] = pwelch(recoveredSine, hamming(nfftWelch), nfftWelch/2, nfftWelch, Fs, 'centered');
plot(fRec/1e3, 10*log10(pxxRec + eps));
grid on;
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('Recovered sine spectrum');

save('SINsignal_recovered_analysis.mat', ...
    'Fs', 'toneFreq', 'complexAmp', 'recoveredSine', 'demodEnvelope', '-v7.3');
fprintf('  Saved recovered result: SINsignal_recovered_analysis.mat\n');

function matFile = locate_mat_file(fileName)
if exist(fileName, 'file')
    matFile = fileName;
    return;
end

candidate = fullfile(pwd, '0622', fileName);
if exist(candidate, 'file')
    matFile = candidate;
    return;
end

listing = dir(fullfile(pwd, '**', fileName));
if isempty(listing)
    error('Cannot find %s under %s.', fileName, pwd);
end
matFile = fullfile(listing(1).folder, listing(1).name);
end

function x = load_first_signal(matFile)
s = load(matFile);
names = fieldnames(s);
v = s.(names{1});

if isa(v, 'timeseries')
    x = v.Data;
elseif isstruct(v) && isfield(v, 'signals') && isfield(v.signals, 'values')
    x = v.signals.values;
elseif isstruct(v) && isfield(v, 'Data')
    x = v.Data;
else
    x = v;
end
end

function [toneFreq, tonePowerDb] = estimate_tone_frequency(x, Fs)
nfft = 2^nextpow2(min(numel(x), 1048576));
seg = x(1:min(numel(x), nfft));
win = hann(numel(seg));
spec = fftshift(fft(seg(:).*win, nfft));
freq = (-nfft/2:nfft/2-1).' * Fs/nfft;
powerDb = 20*log10(abs(spec) + eps);

dcGuard = abs(freq) < 100;
powerDb(dcGuard) = -Inf;
[tonePowerDb, idx] = max(powerDb);
toneFreq = freq(idx);

if idx > 1 && idx < numel(freq)
    y1 = powerDb(idx-1);
    y2 = powerDb(idx);
    y3 = powerDb(idx+1);
    denom = y1 - 2*y2 + y3;
    if isfinite(denom) && abs(denom) > eps
        delta = 0.5 * (y1 - y3) / denom;
        toneFreq = toneFreq + delta * Fs/nfft;
    end
end
end
