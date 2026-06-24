function rxSym = equalize_payload(dataIn, pilotIn)
%#codegen

M = 4;
numDataSC = 48;
numPayloadSym = 10;

rng(2);
trainIdx = randi([0 M-1], numDataSC, 1);
trainSym = pskmod(trainIdx, M, pi/4, 'gray');

H = dataIn(:,1) ./ trainSym;
Y = dataIn(:,2:end) ./ H;

Hpilot = pilotIn(:,1);

for k = 1:numPayloadSym
    cpe = mean(pilotIn(:,k+1) ./ Hpilot);
    if abs(cpe) > 0
        Y(:,k) = Y(:,k) ./ (cpe / abs(cpe));
    end
end

rxSym = double(Y(:));