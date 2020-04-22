function maxAmp = rescaleFigWav(hFigWav, hClust, maxAmpPrev, scaleFactor, channel_idx)
    %RESCALEFIGWAV
    maxAmp = maxAmpPrev * scaleFactor;
    jrclust.views.plotFigWav(hFigWav, hClust, maxAmp, channel_idx);
end
