function updateCursorFigSim(obj)
    %UPDATECURSORFIGSIM Update the crosshair cursor in FigSim
    if isempty(obj.selected) || ~obj.hasFig('FigSim')
        return;
    end

    iCluster = obj.selected(1);
    if numel(obj.selected) == 1
        jCluster = iCluster;
    else
        jCluster = obj.selected(2);
    end

    hFigSim = obj.hFigs('FigSim');

    % update crosshair cursor
    hFigSim.updatePlot('hCursorV', iCluster*[1, 1], 0.5 + [0, obj.hClust.nClusters]);
    if iCluster == jCluster
        colorH = obj.hCfg.colorMap(2, :); % black
    else
        colorH = obj.hCfg.colorMap(3, :); % red
    end
    hFigSim.updatePlot('hCursorH', 0.5 + [0, obj.hClust.nClusters], jCluster*[1, 1]);
    hFigSim.plotApply('hCursorH', @set, 'Color', colorH);

    % center on this pair of clusters
    hFigSim.axApply('default', @set, 'XLim', jrclust.utils.trimLim(iCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));
    hFigSim.axApply('default', @set, 'YLim', jrclust.utils.trimLim(jCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));

    if strcmp(hFigSim.figData.figView, 'template')
        scoreij = obj.hClust.templateSim(iCluster, jCluster);
    elseif strcmp(hFigSim.figData.figView, 'waveform')
        scoreij = obj.hClust.waveformSim(iCluster, jCluster);
    end
    hFigSim.axApply('default', @title, sprintf('Cluster %d vs. Cluster %d: %0.3f', iCluster, jCluster, scoreij));
end