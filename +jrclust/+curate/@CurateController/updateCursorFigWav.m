function updateCursorFigWav(obj)
    %UPDATECURSORFIGWAV Plot mean waveforms on top of selected cluster(s)
    if isempty(obj.selected) || ~obj.hasFig('FigWav')
        return;
    end
    hFigWav = obj.hFigs('FigWav');
    plotSelectedMeansFun = @(x,y)plotSelectedMeans(hFigWav, obj.hClust, x, y, obj.maxAmp, obj.hCfg, obj.spatial_idx);
    if numel(obj.selected) == 1 % we've selected just one cluster (primary)
        hFigWav.rmPlot('selected2'); % if already selected, hide it
        plotSelectedMeansFun(obj.selected,'selected1');
    else % just plot #2 for now
        plotSelectedMeansFun(obj.selected(1),'selected1');
        plotSelectedMeansFun(obj.selected(2),'selected2');
    end
end

%% LOCAL FUNCTIONS
function iCluster = plotSelectedMeans(hFigWav, hClust, iCluster, plotKey, maxAmp, hCfg, spatial_idx)
    %PLOTSELECTEDMEANS Plot an overlay on selected cluster's mean traces
    showSubset = hFigWav.figData.showSubset;

    if strcmp(plotKey, 'selected2')
        colorMap = hCfg.colorMap(3, :); % red
    elseif strcmp(plotKey, 'selected1')
        colorMap = hCfg.colorMap(2, :); % black
    else
        return; % maybe some day we support selected3, 4, ...
    end

    if ~hFigWav.hasPlot(plotKey)
        hFigWav.addPlot(plotKey, nan, nan, 'Color', colorMap, 'LineWidth', 2);
    end

    if hCfg.showRaw
        meanWf = hClust.meanWfGlobalRaw(:, :, iCluster);
    else
        meanWf = hClust.meanWfGlobal(:, :, iCluster);
    end

    hFigWav.multiplot(plotKey, maxAmp, jrclust.views.getXRange(find(iCluster == showSubset), size(meanWf, 1), hCfg), meanWf(:,spatial_idx));
    hFigWav.plotApply(plotKey, @uistack, 'top');
end
