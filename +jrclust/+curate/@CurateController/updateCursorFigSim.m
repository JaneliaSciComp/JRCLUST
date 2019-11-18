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
    
    iClusterX = obj.unitIndex(iCluster);
    jClusterY = obj.unitIndex(jCluster);

    hFigSim = obj.hFigs('FigSim');

    % update crosshair cursor
    hFigSim.updatePlot('hCursorV', iClusterX*[1, 1], 0.5 + [0, obj.nShown]);
    if iCluster == jCluster
        colorH = obj.hCfg.colorMap(2, :); % black
    else
        colorH = obj.hCfg.colorMap(3, :); % red
    end
    hFigSim.updatePlot('hCursorH', 0.5 + [0, obj.nShown], jClusterY*[1, 1]);
    hFigSim.plotApply('hCursorH', @set, 'Color', colorH);

    % center on this pair of clusters
    hFigSim.axApply('default', @set, 'XLim', jrclust.utils.trimLim(iClusterX + [-6, 6], 0.5 + [0, obj.nShown]));
    hFigSim.axApply('default', @set, 'YLim', jrclust.utils.trimLim(jClusterY + [-6, 6], 0.5 + [0, obj.nShown]));

    if strcmp(hFigSim.figData.figView, 'template')
        scoreij = obj.hClust.templateSim(iCluster, jCluster);
    elseif strcmp(hFigSim.figData.figView, 'waveform')
        scoreij = obj.hClust.waveformSim(iCluster, jCluster);
    end

    if iCluster ~= jCluster
        hFigSim.axApply('default', @title, sprintf('Unit %d vs. Unit %d: %0.3f (%s) (press [H] for help)', iCluster, jCluster, scoreij, hFigSim.figData.figView));
    else
        hFigSim.axApply('default', @title, sprintf('Unit %d: %0.3f (%s) (press [H] for help)', iCluster, scoreij, hFigSim.figData.figView));
    end
end