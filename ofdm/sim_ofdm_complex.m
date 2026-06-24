% Simulation-level OFDM link matching the original Simulink plan.
% Features: preamble sync, coarse CFO correction, pilot-aided CPE correction,
% one training OFDM symbol, ten payload OFDM symbols, QPSK, BER measurement.

clear; clc; close all;

cfg.Fs = 1e6;
cfg.Nfft = 64;
cfg.cpLen = 16;
cfg.guardBands = [6; 5];
cfg.pilotIdx = [12; 26; 40; 54];
cfg.M = 4;
cfg.phaseOffset = pi/4;
cfg.numTrainSym = 1;
cfg.numPayloadSym = 10;
cfg.numOFDMSym = cfg.numTrainSym + cfg.numPayloadSym;
cfg.txScale = 0.15;
cfg.snrDb = 30;
cfg.cfoHz = 1200;
cfg.syncThreshold = 0.45;

idx = carrier_indices(cfg.Nfft, cfg.guardBands, cfg.pilotIdx);
numDataSC = numel(idx.data);

rng(1);
preHalf = exp(1j*2*pi*(0:63).'/8);
preamble = [preHalf; preHalf];

rng(2);
trainIdx = randi([0 cfg.M-1], numDataSC, 1);
payloadIdx = randi([0 cfg.M-1], numDataSC*cfg.numPayloadSym, 1);

trainSym = pskmod(trainIdx, cfg.M, cfg.phaseOffset, 'gray');
payloadSym = reshape(pskmod(payloadIdx, cfg.M, cfg.phaseOffset, 'gray'), ...
    numDataSC, cfg.numPayloadSym);
pilotSym = ones(numel(cfg.pilotIdx), cfg.numOFDMSym);

txDataGrid = [trainSym, payloadSym];
ofdmTx = ofdm_modulate(txDataGrid, pilotSym, cfg, idx);
txFrame = cfg.txScale * [preamble; ofdmTx];

rx = impair_channel(txFrame, cfg);
[ofdmFrame, valid, syncMetric, startIdx, estCfoHz] = sync_and_correct(rx, preamble, cfg);
if ~valid
    error('Frame sync failed. Best metric %.3f is below threshold %.3f.', ...
        syncMetric, cfg.syncThreshold);
end

[dataOut, pilotOut] = ofdm_demodulate(ofdmFrame, cfg, idx);

Hdata = dataOut(:,1) ./ trainSym;
payloadEq = dataOut(:,2:end) ./ Hdata;

Hpilot = pilotOut(:,1) ./ pilotSym(:,1);
for k = 1:cfg.numPayloadSym
    cpe = mean((pilotOut(:,k+1) ./ pilotSym(:,k+1)) ./ Hpilot);
    if abs(cpe) > 0
        payloadEq(:,k) = payloadEq(:,k) ./ (cpe / abs(cpe));
    end
end

rxPayloadSym = payloadEq(:);
rxIdx = pskdemod(rxPayloadSym, cfg.M, cfg.phaseOffset, 'gray');

numErr = nnz(rxIdx ~= payloadIdx);
ber = numErr / numel(payloadIdx);

fprintf('Complex OFDM simulation\n');
fprintf('  Data subcarriers: %d\n', numDataSC);
fprintf('  Frame sync metric: %.3f\n', syncMetric);
fprintf('  Frame start index: %d\n', startIdx);
fprintf('  True CFO: %.1f Hz, estimated CFO: %.1f Hz\n', cfg.cfoHz, estCfoHz);
fprintf('  BER: %.6g (%d / %d symbol errors)\n', ber, numErr, numel(payloadIdx));

figure('Name', 'Complex OFDM Equalized QPSK');
plot(real(rxPayloadSym), imag(rxPayloadSym), '.');
grid on; axis equal;
xlabel('In-phase'); ylabel('Quadrature');
title(sprintf('Complex OFDM, BER = %.3g', ber));

function idx = carrier_indices(Nfft, guardBands, pilotIdx)
dcIdx = Nfft/2 + 1;
guardIdx = [1:guardBands(1), Nfft-guardBands(2)+1:Nfft].';
nullIdx = unique([guardIdx; dcIdx; pilotIdx(:)]);
dataIdx = setdiff((1:Nfft).', nullIdx, 'stable');

idx.dc = dcIdx;
idx.guard = guardIdx;
idx.pilot = pilotIdx(:);
idx.data = dataIdx;
end

function y = ofdm_modulate(dataGrid, pilotGrid, cfg, idx)
numSym = size(dataGrid, 2);
y = complex(zeros(numSym*(cfg.Nfft+cfg.cpLen), 1));

for k = 1:numSym
    freqCentered = complex(zeros(cfg.Nfft, 1));
    freqCentered(idx.data) = dataGrid(:,k);
    freqCentered(idx.pilot) = pilotGrid(:,k);

    timeSym = ifft(ifftshift(freqCentered)) * sqrt(cfg.Nfft);
    withCp = [timeSym(end-cfg.cpLen+1:end); timeSym];

    outIdx = (k-1)*(cfg.Nfft+cfg.cpLen) + (1:(cfg.Nfft+cfg.cpLen));
    y(outIdx) = withCp;
end
end

function [dataOut, pilotOut] = ofdm_demodulate(x, cfg, idx)
dataOut = complex(zeros(numel(idx.data), cfg.numOFDMSym));
pilotOut = complex(zeros(numel(idx.pilot), cfg.numOFDMSym));

symLen = cfg.Nfft + cfg.cpLen;
for k = 1:cfg.numOFDMSym
    inIdx = (k-1)*symLen + cfg.cpLen + (1:cfg.Nfft);
    freqCentered = fftshift(fft(x(inIdx)) / sqrt(cfg.Nfft));
    dataOut(:,k) = freqCentered(idx.data);
    pilotOut(:,k) = freqCentered(idx.pilot);
end
end

function rx = impair_channel(tx, cfg)
n = (0:numel(tx)-1).';
cfoNorm = cfg.cfoHz / cfg.Fs;

channelTaps = [1; 0.20*exp(1j*0.65); 0.08*exp(-1j*1.1)];
rx = filter(channelTaps, 1, tx);
rx = rx .* exp(1j*2*pi*cfoNorm*n);
rx = add_awgn(rx, cfg.snrDb);
end

function y = add_awgn(x, snrDb)
sigPower = mean(abs(x).^2);
noisePower = sigPower / 10^(snrDb/10);
noise = sqrt(noisePower/2) * (randn(size(x)) + 1j*randn(size(x)));
y = x + noise;
end

function [ofdmFrame, valid, bestMetric, bestIdx, estCfoHz] = sync_and_correct(rx, preamble, cfg)
ofdmLen = cfg.numOFDMSym * (cfg.Nfft + cfg.cpLen);
ofdmFrame = complex(zeros(ofdmLen, 1));
valid = false;
bestMetric = 0;
bestIdx = 1;
estCfoHz = 0;

rx = rx(:);
searchLast = numel(rx) - numel(preamble) - ofdmLen + 1;
if searchLast < 1
    return;
end

preNorm = sqrt(sum(abs(preamble).^2));
for k = 1:searchLast
    seg = rx(k:k+numel(preamble)-1);
    metric = abs(sum(conj(preamble).*seg)) / ...
        (sqrt(sum(abs(seg).^2))*preNorm + eps);
    if metric > bestMetric
        bestMetric = metric;
        bestIdx = k;
    end
end

if bestMetric < cfg.syncThreshold
    return;
end

r1 = rx(bestIdx:bestIdx+63);
r2 = rx(bestIdx+64:bestIdx+127);
phi = angle(sum(conj(r1).*r2));
cfoNorm = phi / (2*pi*cfg.Nfft);
estCfoHz = cfoNorm * cfg.Fs;

startIdx = bestIdx + numel(preamble);
x = rx(startIdx:startIdx+ofdmLen-1);
n = (0:ofdmLen-1).';
ofdmFrame = x .* exp(-1j*2*pi*cfoNorm*n);
valid = true;
end
