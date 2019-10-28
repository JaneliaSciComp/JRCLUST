function maxAmp = rescaleFigWav(hFigWav, hClust, maxAmpPrev, scaleFactor)
    %RESCALEFIGWAV
    maxAmp = maxAmpPrev * scaleFactor;
    jrclust.views.plotFigWav(hFigWav, hClust, maxAmp);
end
