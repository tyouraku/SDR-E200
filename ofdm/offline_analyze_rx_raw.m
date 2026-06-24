% Offline analysis for saved Pluto/AntSDR raw IQ data.
% This script reads ofdm_rx_raw.mat without modifying any existing files.

clear; clc; close all;

matFile = locate_mat_file('ofdm_rx_raw.mat');
Fs = 1e6;
Nfft = 64;
cpLen = 16;
numOFDMSym = 11;
ofdmLen = numOFDMSym * (Nfft + cpLen);

rxCube = load_first_signal(matFile);
rxRaw = rxCube(:);
rxRaw = rxRaw(~isnan(real(rxRaw)) & ~isnan(imag(rxRaw)));

rng(1);
preHalf = exp(1j*2*pi*(0:63).'/8);
preamble = [preHalf; preHalf];

[metric, peakIdx, peakVal] = preamble_metric(rxRaw, preamble, ofdmLen);
estCfoHz = estimate_cfo(rxRaw, peakIdx, Fs, Nfft);

fprintf('Raw IQ offline analysis\n');
fprintf('  File: %s\n', matFile);
fprintf('  Total samples: %d\n', numel(rxRaw));
fprintf('  RMS amplitude: %.4g\n', rms(abs(rxRaw)));
fprintf('  Peak amplitude: %.4g\n', max(abs(rxRaw)));
fprintf('  PAPR: %.2f dB\n', 20*log10(max(abs(rxRaw))/(rms(abs(rxRaw)) + eps)));
fprintf('  Best preamble metric: %.3f at sample %d\n', peakVal, peakIdx);
fprintf('  Estimated CFO from preamble: %.1f Hz\n', estCfoHz);

numPlot = min(numel(rxRaw), 20000);
t = (0:numPlot-1).' / Fs;

figure('Name', 'Raw IQ Time Domain');
subplot(2,1,1);
plot(t, real(rxRaw(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('I');
title('Raw IQ real part');

subplot(2,1,2);
plot(t, abs(rxRaw(1:numPlot)));
grid on;
xlabel('Time (s)'); ylabel('|r[n]|');
title('Raw IQ amplitude');

figure('Name', 'Raw IQ Spectrum');
nfftWelch = 4096;
[pxx, f] = pwelch(rxRaw, hamming(nfftWelch), nfftWelch/2, nfftWelch, Fs, 'centered');
plot(f/1e3, 10*log10(pxx + eps));
grid on;
xlabel('Frequency (kHz)'); ylabel('PSD (dB/Hz)');
title('Raw IQ spectrum');

figure('Name', 'Preamble Correlation Metric');
plot(metric);
grid on;
xlabel('Sample index'); ylabel('Normalized correlation');
title(sprintf('Preamble metric, peak %.3f at %d', peakVal, peakIdx));

if peakIdx + numel(preamble) + ofdmLen - 1 <= numel(rxRaw)
    oneFrame = rxRaw(peakIdx:peakIdx+numel(preamble)+ofdmLen-1);
    figure('Name', 'Detected Frame Amplitude');
    plot(abs(oneFrame));
    grid on;
    xlabel('Sample in detected frame'); ylabel('Amplitude');
    title('Detected preamble + OFDM frame amplitude');
end

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

function [metric, peakIdx, peakVal] = preamble_metric(rx, preamble, ofdmLen)
searchLast = numel(rx) - numel(preamble) - ofdmLen + 1;
if searchLast < 1
    error('Signal is shorter than one preamble + OFDM frame.');
end

metric = zeros(searchLast, 1);
preNorm = sqrt(sum(abs(preamble).^2));
for k = 1:searchLast
    seg = rx(k:k+numel(preamble)-1);
    metric(k) = abs(sum(conj(preamble).*seg)) / ...
        (sqrt(sum(abs(seg).^2))*preNorm + eps);
end

[peakVal, peakIdx] = max(metric);
end

function estCfoHz = estimate_cfo(rx, peakIdx, Fs, Nfft)
if peakIdx + 127 > numel(rx)
    estCfoHz = NaN;
    return;
end
r1 = rx(peakIdx:peakIdx+63);
r2 = rx(peakIdx+64:peakIdx+127);
phi = angle(sum(conj(r1).*r2));
estCfoHz = phi / (2*pi*Nfft) * Fs;
end
