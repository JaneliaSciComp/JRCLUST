function maxAmp = rescaleFigWav(hFigWav, hClust, hCfg, maxAmpPrev, scaleFactor)
    %RESCALEFIGWAV
    maxAmp = maxAmpPrev * scaleFactor;
    doPlotFigWav(hFigWav, hClust, hCfg, maxAmp);
end
