function updateFigWav(obj)
    %UPDATEFIGWAV
    if ~obj.hasFig('FigWav')
        return;
    end
    fprintf('Number of clusters at update FigWav: %d\n', obj.hClust.nClusters);
    
    hFigWav = obj.hFigs('FigWav');
    obj.showSubset = 1:obj.hClust.nClusters;
    jrclust.views.plotFigWav(hFigWav, obj.hClust, obj.maxAmp, obj.showSubset, obj.channel_idx);
    setFigWavXTicks(hFigWav, obj.hClust, obj.hCfg.showSpikeCount);
end
