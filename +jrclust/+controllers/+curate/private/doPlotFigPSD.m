%function [vrPow, vrFreq] = doPlotFigPSD(hFigPSD, tracesFilt, hCfg)
function hFigPSD = doPlotFigPSD(hFigPSD, tracesFilt, hCfg)
    nSmooth = 3;
    ignoreSites = hCfg.ignoreSites;
    sampleRate = hCfg.sampleRate/hCfg.nSkip_show;

    %warning off;

    tracesFilt = fft(tracesFilt');
    tracesFilt = real(tracesFilt .* conj(tracesFilt)) / size(tracesFilt,1) / (sampleRate/2);

    imid0 = ceil(size(tracesFilt, 1)/2);
    vrFreq = (0:size(tracesFilt, 1) - 1)*(sampleRate/size(tracesFilt, 1));
    vrFreq = vrFreq(2:imid0);
    viChan = setdiff(1:size(tracesFilt,2), ignoreSites);

    if size(tracesFilt, 2) > 1
        vrPow = mean(tracesFilt(2:imid0,viChan), 2);
        YLabel = 'Mean power across sites (dB uV^2/Hz)';
    else
        vrPow = tracesFilt(2:imid0, 1);
        YLabel = 'Power on site (dB uV^2/Hz)';
    end

    vrPow = filterq_(ones([nSmooth, 1]), nSmooth, vrPow);

    hFigPSD.addPlot('hFreq', vrFreq, jrclust.utils.pow2db(vrPow), 'k-');
    %     set(gca, 'YScale', 'log');
    %     set(gca, 'XScale', 'linear');
    hFigPSD.axApply('default', @xlabel, 'Frequency (Hz)');
    hFigPSD.axApply('default', @ylabel, YLabel);
    % xlim_([0 sRateHz/2]);
    hFigPSD.axApply('default', @grid, 'on');

    try
        hFigPSD.axApply('default', @xlim, vrFreq([1, end]));
        hFigPSD.figApply(@set, 'Color', 'w');
    catch
    end
end
