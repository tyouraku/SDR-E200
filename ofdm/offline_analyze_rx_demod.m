% Offline analysis for saved OFDM demodulator output.
% This script reads ofdm_rx_demod.mat without modifying any existing files.

clear; clc; close all;

matFile = locate_mat_file('ofdm_rx_demod.mat');
M = 4;
phaseOffset = pi/4;

dataCube = load_first_signal(matFile);
dataCube = squeeze(dataCube);

if ndims(dataCube) == 2
    dataCube = reshape(dataCube, size(dataCube,1), size(dataCube,2), 1);
end

[numDataSC, numOFDMSym, numFrames] = size(dataCube);

rng(2);
trainIdx = randi([0 M-1], numDataSC, 1);
trainSym = pskmod(trainIdx, M, phaseOffset, 'gray');

payloadIdx = [];
if numOFDMSym > 1
    payloadIdx = randi([0 M-1], numDataSC*(numOFDMSym-1), 1);
end

H = complex(zeros(numDataSC, numFrames));
rxEq = complex(zeros(numDataSC, max(numOFDMSym-1, 0), numFrames));
serPerFrame = nan(numFrames, 1);
evmPerFrame = nan(numFrames, 1);
validFrame = false(numFrames, 1);

for f = 1:numFrames
    frameDemod = dataCube(:,:,f);
    frameFinite = isfinite(real(frameDemod)) & isfinite(imag(frameDemod));
    if ~all(frameFinite(:))
        continue;
    end

    H(:,f) = frameDemod(:,1) ./ trainSym;
    goodH = isfinite(real(H(:,f))) & isfinite(imag(H(:,f))) & abs(H(:,f)) > 1e-8;
    if nnz(goodH) < 0.8*numDataSC
        continue;
    end

    if numOFDMSym > 1
        eqFrame = frameDemod(:,2:end) ./ H(:,f);
        rxEq(:,:,f) = eqFrame;

        eqVec = eqFrame(:);
        goodEq = isfinite(real(eqVec)) & isfinite(imag(eqVec));
        if ~any(goodEq)
            continue;
        end

        rxIdx = pskdemod(eqVec(goodEq), M, phaseOffset, 'gray');
        serPerFrame(f) = nnz(rxIdx ~= payloadIdx(goodEq)) / nnz(goodEq);

        refSym = pskmod(payloadIdx, M, phaseOffset, 'gray');
        evmPerFrame(f) = sqrt(mean(abs(eqVec(goodEq) - refSym(goodEq)).^2) / ...
            (mean(abs(refSym(goodEq)).^2) + eps));
        validFrame(f) = true;
    end
end

fprintf('OFDM demod offline analysis\n');
fprintf('  File: %s\n', matFile);
fprintf('  Data cube size: %d subcarriers x %d OFDM symbols x %d frames\n', ...
    numDataSC, numOFDMSym, numFrames);
goodHAll = H(isfinite(real(H)) & isfinite(imag(H)) & abs(H) > 0);
fprintf('  Valid analyzed frames: %d / %d\n', nnz(validFrame), numFrames);
fprintf('  Median |H|: %.4g\n', median(abs(goodHAll)));
fprintf('  Median channel ripple: %.2f dB\n', median(channel_ripple_db(H(:,validFrame)), 'omitnan'));
if numOFDMSym > 1
    fprintf('  Median SER: %.6g\n', median(serPerFrame, 'omitnan'));
    fprintf('  Median EVM: %.2f dB\n', 20*log10(median(evmPerFrame, 'omitnan') + eps));
end

firstGood = find(validFrame, 1, 'first');
if isempty(firstGood)
    firstGood = 1;
end
firstFrame = dataCube(:,:,firstGood);
figure('Name', 'Demodulated Subcarriers Before Equalization');
plot(real(firstFrame(:)), imag(firstFrame(:)), '.');
grid on; axis equal;
xlabel('In-phase'); ylabel('Quadrature');
title(sprintf('OFDM demod output before equalization, frame %d', firstGood));

figure('Name', 'Estimated Channel Response');
subplot(2,1,1);
plot(20*log10(abs(H(:,firstGood)) + eps), 'o-');
grid on;
xlabel('Data subcarrier index'); ylabel('|H| (dB)');
title(sprintf('Estimated channel magnitude, frame %d', firstGood));

subplot(2,1,2);
plot(unwrap(angle(H(:,firstGood))), 'o-');
grid on;
xlabel('Data subcarrier index'); ylabel('Phase (rad)');
title(sprintf('Estimated channel phase, frame %d', firstGood));

if numOFDMSym > 1
    eqAll = rxEq(:);
    eqAll = eqAll(~isnan(real(eqAll)) & ~isnan(imag(eqAll)));

    figure('Name', 'Equalized Payload Constellation');
    plot(real(eqAll), imag(eqAll), '.');
    grid on; axis equal;
    xlabel('In-phase'); ylabel('Quadrature');
    title('Equalized QPSK payload constellation');

    figure('Name', 'Frame Quality Metrics');
    subplot(2,1,1);
    plot(serPerFrame, 'o-');
    grid on;
    xlabel('Frame index'); ylabel('SER');
    title('Symbol error rate per saved frame');

    subplot(2,1,2);
    plot(20*log10(evmPerFrame + eps), 'o-');
    grid on;
    xlabel('Frame index'); ylabel('EVM (dB)');
    title('EVM per saved frame');
end

figure('Name', 'Channel Magnitude Over Time');
imagesc(20*log10(abs(H) + eps));
axis xy; colorbar;
xlabel('Frame index'); ylabel('Data subcarrier index');
title('Estimated channel magnitude over saved frames');

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

function r = channel_ripple_db(H)
r = zeros(size(H,2), 1);
for k = 1:size(H,2)
    magDb = 20*log10(abs(H(:,k)) + eps);
    r(k) = max(magDb) - min(magDb);
end
end
