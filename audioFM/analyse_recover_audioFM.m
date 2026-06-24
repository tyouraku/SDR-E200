%% redo_recover_audioFM.m
% Fresh recovery flow for D:\drone\0511\audioFM.mat.
% This script does not modify the original project files.

clc;
clear;
close all;

cfg = struct();
cfg.workDir = 'D:\drone\0511\重做恢复';
cfg.inMat = fullfile(cfg.workDir, 'audioFMnew1.mat');
cfg.refAudio = 'D:\drone\Matlab\audio_sources\music1_mono48kHz.wav';
cfg.outDir = fullfile(cfg.workDir, 'out');

% Parameters copied from the Simulink chain.
cfg.rxFs = 240e3;
cfg.fmDeviationHz = 75e3;
cfg.rfAudioPassbandHz = 15e3;
cfg.rfAudioStopbandHz = 20e3;
cfg.audioFs = 48e3;
cfg.decimation = 5;

% Recovery parameters for the new attempt.
cfg.hpCutoffHz = 60;
cfg.enableDeemphasis = false;
cfg.deemphasisTauSec = 75e-6;
cfg.enableHumNotch = true;
cfg.humCandidatesHz = [50 60];
cfg.maxHumHarmonic = 4;
cfg.enableSpectralDenoise = false;
cfg.stftLength = 1024;
cfg.stftOverlap = 768;
cfg.stftFloor = 0.12;
cfg.noisePercentile = 0.15;
cfg.overSubtract = 1.10;
cfg.targetPeak = 0.97;
cfg.normPercentile = 99.5;
cfg.maxNormGain = 4.0;
cfg.analysisSpanSec = 0.08;
cfg.showFigures = true;

if ~isfolder(cfg.outDir)
    mkdir(cfg.outDir);
end

summaryPath = fullfile(cfg.outDir, 'recovery_summary.txt');
fid = fopen(summaryPath, 'w');
cleanupObj = onCleanup(@() fclose(fid));

logLine(fid, 'Fresh recovery run: %s', datestr(now, 31));
logLine(fid, 'Input MAT: %s', cfg.inMat);
logLine(fid, 'Output folder: %s', cfg.outDir);

raw = load(cfg.inMat);
[iq, fsRx, meta] = loadRxIq(raw, cfg.rxFs);
logLine(fid, 'Loaded variable: %s', meta.varName);
logLine(fid, 'Detected fsRx: %.3f Hz', fsRx);
logLine(fid, 'Signal shape: %s', meta.shapeText);
logLine(fid, 'Frame length: %d', meta.frameLen);
logLine(fid, 'Frame count: %d', meta.frameCount);
logLine(fid, 'Duration: %.3f s', numel(iq) / fsRx);

rawStats = analyzeRealImag(iq, fsRx);
writeStats(fid, 'RX IQ stats', rawStats);
plotIqDiagnostics(iq, fsRx, cfg, fullfile(cfg.outDir, '01_rx_iq_analysis.png'));

[audioDemod240k, instFreqHz, freqOffsetHz] = fmDemodFresh(iq, fsRx, cfg.fmDeviationHz);
logLine(fid, 'Estimated center offset from demod median: %.3f Hz', freqOffsetHz);

plotDemodDiagnostics(instFreqHz, audioDemod240k, fsRx, cfg, ...
    fullfile(cfg.outDir, '02_demod_240k_analysis.png'));

audioDemod48k = postDemodChain(audioDemod240k, fsRx, cfg);
demodStats = analyzeAudio(audioDemod48k, cfg.audioFs);
writeStats(fid, 'Recovered audio before enhancement', demodStats);
plotAudioDiagnostics(audioDemod48k, cfg.audioFs, cfg, ...
    'Recovered Audio Before Enhancement', ...
    fullfile(cfg.outDir, '03_audio_before_enhance.png'));

[audioOpt48k, enhanceInfo] = enhanceRecoveredAudio(audioDemod48k, cfg);
optStats = analyzeAudio(audioOpt48k, cfg.audioFs);
writeStats(fid, 'Recovered audio after enhancement', optStats);
logLine(fid, 'Detected hum fundamental: %.1f Hz', enhanceInfo.humFundamentalHz);
logLine(fid, 'Applied notch frequencies: %s', mat2str(enhanceInfo.notchFreqsHz));

plotAudioDiagnostics(audioOpt48k, cfg.audioFs, cfg, ...
    'Recovered Audio After Enhancement', ...
    fullfile(cfg.outDir, '04_audio_after_enhance.png'));
plotBeforeAfter(audioDemod48k, audioOpt48k, cfg.audioFs, cfg, ...
    fullfile(cfg.outDir, '05_before_after_compare.png'));

refMetrics = struct();
if isfile(cfg.refAudio)
    [refMetrics, refWave] = compareWithReference(audioOpt48k, cfg.audioFs, cfg.refAudio); %#ok<ASGLU>
    writeReferenceMetrics(fid, refMetrics);
    plotReferenceCompare(refWave.refAligned, refWave.recAligned, cfg.audioFs, cfg, ...
        fullfile(cfg.outDir, '06_reference_compare.png'));
else
    logLine(fid, 'Reference audio not found: %s', cfg.refAudio);
end

demodWav = fullfile(cfg.outDir, 'audio_demod_48k.wav');
optWav = fullfile(cfg.outDir, 'audio_recovered_optimized_48k.wav');
audiowrite(demodWav, audioDemod48k, cfg.audioFs);
audiowrite(optWav, audioOpt48k, cfg.audioFs);

resultPath = fullfile(cfg.outDir, 'redo_recovery_resultnew1.mat');
save(resultPath, 'cfg', 'meta', 'iq', 'fsRx', 'rawStats', 'instFreqHz', ...
    'audioDemod240k', 'audioDemod48k', 'audioOpt48k', 'demodStats', ...
    'optStats', 'enhanceInfo', 'refMetrics', 'freqOffsetHz', '-v7.3');

logLine(fid, 'Saved demod wav: %s', demodWav);
logLine(fid, 'Saved optimized wav: %s', optWav);
logLine(fid, 'Saved result mat: %s', resultPath);
logLine(fid, 'Done.');

fprintf('Fresh recovery completed.\n');
fprintf('Summary: %s\n', summaryPath);
fprintf('Demod wav: %s\n', demodWav);
fprintf('Optimized wav: %s\n', optWav);

function [iq, fs, meta] = loadRxIq(S, fsFallback)
    names = fieldnames(S);
    picked = '';
    iq = [];
    fs = [];
    meta = struct('varName', '', 'shapeText', '', 'frameLen', 0, 'frameCount', 0);

    for k = 1:numel(names)
        v = S.(names{k});
        if isa(v, 'timeseries')
            data = v.Data;
            sz = size(data);
            if ndims(data) == 3 && sz(2) == 1
                frameLen = sz(1);
                frameCount = sz(3);
                iq = reshape(data(:, 1, :), frameLen * frameCount, 1);
                meta.shapeText = sprintf('timeseries [%d x 1 x %d]', frameLen, frameCount);
                meta.frameLen = frameLen;
                meta.frameCount = frameCount;
                if numel(v.Time) > 1
                    dt = median(diff(v.Time), 'omitnan');
                    if isfinite(dt) && dt > 0
                        fs = frameLen / dt;
                    end
                end
                picked = names{k};
                break;
            elseif isvector(data)
                iq = data(:);
                meta.shapeText = sprintf('timeseries vector [%d x 1]', numel(iq));
                meta.frameLen = numel(iq);
                meta.frameCount = 1;
                if numel(v.Time) > 1
                    dt = median(diff(v.Time), 'omitnan');
                    if isfinite(dt) && dt > 0
                        fs = 1 / dt;
                    end
                end
                picked = names{k};
                break;
            end
        end
    end

    if isempty(iq)
        error('No supported timeseries IQ variable was found in the MAT file.');
    end
    if isempty(fs) || ~isfinite(fs) || fs <= 0
        fs = fsFallback;
    end

    meta.varName = picked;
end

function st = analyzeRealImag(x, fs)
    xr = real(x(:));
    xi = imag(x(:));
    amp = abs(x(:));
    st = struct();
    st.fs = fs;
    st.N = numel(x);
    st.durationSec = numel(x) / fs;
    st.realRms = rms(xr);
    st.imagRms = rms(xi);
    st.ampMean = mean(amp);
    st.ampStd = std(amp);
    st.peak = max(amp);
end

function st = analyzeAudio(x, fs)
    x = x(:);
    st = struct();
    st.fs = fs;
    st.N = numel(x);
    st.durationSec = numel(x) / fs;
    st.rms = rms(x);
    st.peak = max(abs(x));
    smooth = movmean(x, 9);
    residual = x - smooth;
    st.noiseRms = rms(residual);
    st.snrEstimateDb = 20 * log10(max(st.rms, eps) / max(st.noiseRms, eps));
end

function [audioNorm, instHz, freqOffsetHz] = fmDemodFresh(iq, fs, deviationHz)
    iq = iq(:);
    mag = abs(iq);
    iq = iq ./ max(mag, eps);
    phDiff = angle(conj(iq(1:end-1)) .* iq(2:end));
    instHz = [0; phDiff] * fs / (2 * pi);
    freqOffsetHz = median(instHz, 'omitnan');
    instHz = instHz - freqOffsetHz;
    audioNorm = instHz / max(deviationHz, eps);
    audioNorm = audioNorm - mean(audioNorm, 'omitnan');
end

function y48k = postDemodChain(x240k, fsRx, cfg)
    d = designfilt('lowpassfir', ...
        'PassbandFrequency', cfg.rfAudioPassbandHz, ...
        'StopbandFrequency', cfg.rfAudioStopbandHz, ...
        'PassbandRipple', 0.1, ...
        'StopbandAttenuation', 80, ...
        'SampleRate', fsRx);
    y = filtfilt(d, x240k(:));
    [p, q] = rat(cfg.audioFs / fsRx, 1e-12);
    y = resample(y, p, q);
    y = y - mean(y, 'omitnan');
    y48k = normalizeToTargetPeak(y, cfg.targetPeak, cfg.normPercentile);
end

function [y, info] = enhanceRecoveredAudio(x, cfg)
    y = x(:);
    info = struct();

    if cfg.hpCutoffHz > 0
        hp = designfilt('highpassiir', ...
            'FilterOrder', 6, ...
            'HalfPowerFrequency', cfg.hpCutoffHz, ...
            'SampleRate', cfg.audioFs);
        y = filtfilt(hp, y);
    end

    if cfg.enableDeemphasis
        y = applyDeemphasis(y, cfg.audioFs, cfg.deemphasisTauSec);
    end

    info.humFundamentalHz = 0;
    info.notchFreqsHz = [];
    if cfg.enableHumNotch
        humHz = detectHumFundamental(y, cfg.audioFs, cfg.humCandidatesHz);
        if humHz > 0
            [y, notchFreqsHz] = applyHumNotches(y, cfg.audioFs, humHz, cfg.maxHumHarmonic);
            info.humFundamentalHz = humHz;
            info.notchFreqsHz = notchFreqsHz;
        end
    end

    if cfg.enableSpectralDenoise
        yCand = spectralSubtractDenoise(y, cfg.audioFs, cfg);
        if rms(yCand) > 0.15 * rms(y) && max(abs(yCand)) > 1e-4
            y = yCand;
        end
    end

    y = normalizeToTargetPeak(y, cfg.targetPeak, cfg.normPercentile);
end

function y = applyDeemphasis(x, fs, tauSec)
    alpha = exp(-1 / (fs * tauSec));
    b = 1 - alpha;
    a = [1 -alpha];
    y = filter(b, a, x(:));
end

function humHz = detectHumFundamental(x, fs, candidatesHz)
    humHz = 0;
    [pxx, f] = pwelch(x, hamming(8192), 4096, 8192, fs, 'onesided');
    localBand = (f >= 40) & (f <= 70);
    if ~any(localBand)
        return;
    end

    score = zeros(size(candidatesHz));
    for k = 1:numel(candidatesHz)
        c = candidatesHz(k);
        band = (f >= c - 1.5) & (f <= c + 1.5);
        side = ((f >= c - 6) & (f <= c - 3)) | ((f >= c + 3) & (f <= c + 6));
        if any(band) && any(side)
            score(k) = mean(pxx(band)) / max(mean(pxx(side)), eps);
        end
    end

    [bestScore, idx] = max(score);
    if bestScore > 2.5
        humHz = candidatesHz(idx);
    end
end

function [y, notchFreqsHz] = applyHumNotches(x, fs, humHz, maxHarmonic)
    y = x(:);
    notchFreqsHz = [];
    for h = 1:maxHarmonic
        f0 = humHz * h;
        if f0 >= fs / 2 - 50
            break;
        end
        wo = f0 / (fs / 2);
        bw = wo / 35;
        [b, a] = iirnotch(wo, bw);
        y = filtfilt(b, a, y);
        notchFreqsHz(end+1) = f0; %#ok<AGROW>
    end
end

function y = spectralSubtractDenoise(x, fs, cfg)
    win = hann(cfg.stftLength, 'periodic');
    [S, f, t] = stft(x, fs, 'Window', win, ...
        'OverlapLength', cfg.stftOverlap, ...
        'FFTLength', cfg.stftLength); %#ok<ASGLU>
    mag = abs(S);
    ph = angle(S);
    noiseProfile = prctile(mag, cfg.noisePercentile * 100, 2);
    noiseMat = noiseProfile * ones(1, size(mag, 2));
    cleanMag = mag - cfg.overSubtract * noiseMat;
    cleanMag = max(cleanMag, cfg.stftFloor * mag);
    Sclean = cleanMag .* exp(1j * ph);
    y = istft(Sclean, fs, 'Window', win, ...
        'OverlapLength', cfg.stftOverlap, ...
        'FFTLength', cfg.stftLength);
    y = y(:);
    if numel(y) > numel(x)
        y = y(1:numel(x));
    elseif numel(y) < numel(x)
        y(end+1:numel(x), 1) = 0;
    end
end

function y = normalizeToTargetPeak(x, targetPeak, pct)
    scale = prctile(abs(x), pct);
    if ~isfinite(scale) || scale <= 0
        y = x;
        return;
    end
    gain = targetPeak / scale;
    gain = min(gain, 4.0);
    y = x * gain;
    y = max(min(y, 0.999), -0.999);
end

function [metrics, waves] = compareWithReference(x, fs, refPath)
    [ref, fsRef] = audioread(refPath);
    if size(ref, 2) > 1
        ref = mean(ref, 2);
    end
    if fsRef ~= fs
        [p, q] = rat(fs / fsRef, 1e-12);
        ref = resample(ref, p, q);
    end

    x = x(:);
    ref = ref(:);
    maxLag = min(round(2 * fs), min(numel(ref), numel(x)) - 1);
    [c, lags] = xcorr(ref, x, maxLag, 'none');
    c = c / max(sqrt(sum(ref .^ 2) * sum(x .^ 2)), eps);
    [bestCorr, idx] = max(c);
    lag = lags(idx);

    if lag >= 0
        refAligned = ref(1 + lag:end);
        recAligned = x(1:min(numel(x), numel(refAligned)));
        refAligned = refAligned(1:numel(recAligned));
    else
        recAligned = x(1 - lag:end);
        refAligned = ref(1:min(numel(ref), numel(recAligned)));
        recAligned = recAligned(1:numel(refAligned));
    end

    gain = (refAligned' * recAligned) / max(recAligned' * recAligned, eps);
    recAligned = recAligned * gain;

    err = refAligned - recAligned;
    metrics = struct();
    metrics.refPath = refPath;
    metrics.bestLagSamples = lag;
    metrics.bestLagSec = lag / fs;
    metrics.peakNormalizedCorrelation = bestCorr;
    metrics.alignedRmse = sqrt(mean(err .^ 2));
    metrics.alignedSnrDb = 20 * log10(rms(refAligned) / max(rms(err), eps));

    waves = struct();
    waves.refAligned = refAligned;
    waves.recAligned = recAligned;
end

function plotIqDiagnostics(iq, fs, cfg, savePath)
    N = numel(iq);
    t = (0:N-1).' / fs;
    span = min(N, round(cfg.analysisSpanSec * fs));
    [pxx, f] = pwelch(iq, hamming(4096), 2048, 4096, fs, 'centered');

    if isfield(cfg, 'showFigures') && cfg.showFigures
        figVisible = 'on';
    else
        figVisible = 'off';
    end

    fig = figure('Visible', figVisible, 'Color', 'w', 'Position', [100 100 1200 900]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t(1:span) * 1e3, real(iq(1:span)), 'b');
    hold on;
    plot(t(1:span) * 1e3, imag(iq(1:span)), 'r');
    grid on;
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('RX IQ Waveform');
    legend('I', 'Q');

    nexttile;
    nScat = min(6000, N);
    idx = round(linspace(1, N, nScat));
    plot(real(iq(idx)), imag(iq(idx)), '.');
    axis equal;
    grid on;
    xlabel('I');
    ylabel('Q');
    title('IQ Scatter');

    nexttile;
    plot(f / 1e3, 10 * log10(pxx + eps));
    grid on;
    xlabel('Frequency (kHz)');
    ylabel('PSD (dB)');
    title('RX IQ Spectrum');

    nexttile;
    spectrogram(iq, hamming(1024), 768, 1024, fs, 'centered', 'yaxis');
    title('RX IQ Spectrogram');

    drawnow;
    exportgraphics(fig, savePath, 'Resolution', 150);
    if strcmp(figVisible, 'off')
        close(fig);
    end
end

function plotDemodDiagnostics(instHz, audioNorm, fs, cfg, savePath)
    N = numel(audioNorm);
    t = (0:N-1).' / fs;
    span = min(N, round(cfg.analysisSpanSec * fs));
    [pxx, f] = pwelch(audioNorm, hamming(4096), 2048, 4096, fs, 'onesided');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 900]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t(1:span) * 1e3, instHz(1:span));
    grid on;
    xlabel('Time (ms)');
    ylabel('Hz');
    title('Instantaneous Frequency');

    nexttile;
    plot(t(1:span) * 1e3, audioNorm(1:span));
    grid on;
    xlabel('Time (ms)');
    ylabel('Normalized audio');
    title('FM Demod Output at 240 kS/s');

    nexttile([1 2]);
    plot(f / 1e3, 10 * log10(pxx + eps));
    grid on;
    xlabel('Frequency (kHz)');
    ylabel('PSD (dB)');
    title('Demod Output Spectrum');

    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
end

function plotAudioDiagnostics(x, fs, cfg, ttl, savePath)
    N = numel(x);
    t = (0:N-1).' / fs;
    span = min(N, round(cfg.analysisSpanSec * fs));
    [pxx, f] = pwelch(x, hamming(4096), 2048, 4096, fs, 'onesided');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 900]);
    tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t(1:span) * 1e3, x(1:span));
    grid on;
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title([ttl ' Waveform']);

    nexttile;
    histogram(x, 120);
    grid on;
    xlabel('Amplitude');
    ylabel('Count');
    title([ttl ' Histogram']);

    nexttile;
    plot(f / 1e3, 10 * log10(pxx + eps));
    grid on;
    xlabel('Frequency (kHz)');
    ylabel('PSD (dB)');
    title([ttl ' Spectrum']);

    nexttile;
    spectrogram(x, hamming(1024), 768, 1024, fs, 'yaxis');
    title([ttl ' Spectrogram']);

    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
end

function plotBeforeAfter(x1, x2, fs, cfg, savePath)
    N = min(numel(x1), numel(x2));
    x1 = x1(1:N);
    x2 = x2(1:N);
    t = (0:N-1).' / fs;
    span = min(N, round(cfg.analysisSpanSec * fs));
    [pxx1, f] = pwelch(x1, hamming(4096), 2048, 4096, fs, 'onesided');
    [pxx2, ~] = pwelch(x2, hamming(4096), 2048, 4096, fs, 'onesided');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 900]);
    tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t(1:span) * 1e3, x1(1:span), 'Color', [0.8 0.2 0.2]);
    hold on;
    plot(t(1:span) * 1e3, x2(1:span), 'Color', [0.1 0.3 0.8]);
    grid on;
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Before vs After Enhancement');
    legend('Before', 'After');

    nexttile;
    plot(f / 1e3, 10 * log10(pxx1 + eps), 'Color', [0.8 0.2 0.2]);
    hold on;
    plot(f / 1e3, 10 * log10(pxx2 + eps), 'Color', [0.1 0.3 0.8]);
    grid on;
    xlabel('Frequency (kHz)');
    ylabel('PSD (dB)');
    title('Before vs After Spectrum');
    legend('Before', 'After');

    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
end

function plotReferenceCompare(refAligned, recAligned, fs, cfg, savePath)
    N = min(numel(refAligned), numel(recAligned));
    refAligned = refAligned(1:N);
    recAligned = recAligned(1:N);
    t = (0:N-1).' / fs;
    span = min(N, round(cfg.analysisSpanSec * fs));

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 700]);
    tiledlayout(2,1, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(t(1:span) * 1e3, refAligned(1:span), 'k');
    hold on;
    plot(t(1:span) * 1e3, recAligned(1:span), 'b');
    grid on;
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title('Reference vs Recovered Waveform');
    legend('Reference', 'Recovered');

    nexttile;
    plot(t(1:span) * 1e3, refAligned(1:span) - recAligned(1:span), 'r');
    grid on;
    xlabel('Time (ms)');
    ylabel('Error');
    title('Reference Alignment Error');

    exportgraphics(fig, savePath, 'Resolution', 150);
    close(fig);
end

function writeStats(fid, heading, st)
    logLine(fid, '%s', heading);
    names = fieldnames(st);
    for k = 1:numel(names)
        v = st.(names{k});
        if isnumeric(v) && isscalar(v)
            logLine(fid, '  %s: %.6g', names{k}, v);
        end
    end
end

function writeReferenceMetrics(fid, metrics)
    logLine(fid, 'Reference comparison');
    logLine(fid, '  refPath: %s', metrics.refPath);
    logLine(fid, '  bestLagSamples: %d', metrics.bestLagSamples);
    logLine(fid, '  bestLagSec: %.6f', metrics.bestLagSec);
    logLine(fid, '  peakNormalizedCorrelation: %.6f', metrics.peakNormalizedCorrelation);
    logLine(fid, '  alignedRmse: %.6g', metrics.alignedRmse);
    logLine(fid, '  alignedSnrDb: %.6f', metrics.alignedSnrDb);
end

function logLine(fid, fmt, varargin)
    txt = sprintf(fmt, varargin{:});
    fprintf('%s\n', txt);
    fprintf(fid, '%s\n', txt);
end
