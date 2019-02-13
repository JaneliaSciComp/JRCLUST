function updateFigWav(obj)
    %UPDATEFIGWAV
    if ~obj.hasFig('FigWav')
        return;
    end

    jrclust.views.plotFigWav(obj.hFigs('FigWav'), obj.hClust, obj.hCfg, obj.maxAmp);
end