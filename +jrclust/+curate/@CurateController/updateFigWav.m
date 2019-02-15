function updateFigWav(obj)
    %UPDATEFIGWAV
    if ~obj.hasFig('FigWav')
        return;
    end

    hFigWav = obj.hFigs('FigWav');
    jrclust.views.plotFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp);
    setFigWavXTicks(hFigWav, obj.hClust, obj.hCfg.showSpikeCount);
end