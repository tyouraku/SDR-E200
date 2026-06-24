Fs = 1e6;
Fc = 2.45e9;

Nfft = 64;
cpLen = 16;
guardBands = [6; 5];
pilotIdx = [12; 26; 40; 54];

M = 4;
numDataSC = 48;
numTrainSym = 1;
numPayloadSym = 10;
numOFDMSym = 11;

ofdmLen = numOFDMSym * (Nfft + cpLen);   % 880
preambleLen = 128;
frameLen = preambleLen + ofdmLen;        % 1008

txScale = 0.15;

rng(1);
preHalf = exp(1j*2*pi*(0:63).'/8);
preamble = [preHalf; preHalf];

rng(2);
trainIdx = randi([0 M-1], numDataSC, 1);
payloadIdx = randi([0 M-1], numDataSC*numPayloadSym, 1);

trainSym = pskmod(trainIdx, M, pi/4, 'gray');
pilotSym = ones(4, numOFDMSym);