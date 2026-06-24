function [payload, valid] = sync_cfo(rx)
%#codegen

Nfft = 64;
cpLen = 16;
numOFDMSym = 11;
ofdmLen = numOFDMSym * (Nfft + cpLen);

payload = complex(zeros(ofdmLen,1,'single'));
valid = false;

preHalf = exp(1j*2*pi*(0:63).'/8);
preamble = [preHalf; preHalf];

rx = rx(:);
L = length(rx);

if L < length(preamble) + ofdmLen
    return;
end

bestMetric = single(0);
bestIdx = 1;

for k = 1:(L - length(preamble) - ofdmLen + 1)
    seg = rx(k:k+127);
    metric = abs(sum(conj(single(preamble)) .* seg)) / ...
        (sqrt(sum(abs(seg).^2)) * sqrt(sum(abs(single(preamble)).^2)) + eps('single'));
    if metric > bestMetric
        bestMetric = metric;
        bestIdx = k;
    end
end

if bestMetric < 0.45
    return;
end

r1 = rx(bestIdx:bestIdx+63);
r2 = rx(bestIdx+64:bestIdx+127);

phi = angle(sum(conj(r1).*r2));
cfoNorm = phi / (2*pi*Nfft);

startIdx = bestIdx + 128;
x = rx(startIdx:startIdx+ofdmLen-1);

n = single((0:ofdmLen-1).');
payload = x .* exp(-1j*2*pi*single(cfoNorm)*n);

valid = true;