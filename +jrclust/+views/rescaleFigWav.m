function maxAmp = rescaleFigWav(hFigWav, hClust, hCfg, maxAmpPrev, scaleFactor)
    %RESCALEFIGWAV
    maxAmp = maxAmpPrev * scaleFactor;
    jrclust.views.plotFigWav(hFigWav, hClust, hCfg, maxAmp);
end
