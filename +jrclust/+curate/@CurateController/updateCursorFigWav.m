function updateCursorFigWav(obj)
    %UPDATECURSORFIGWAV Plot mean waveforms on top of selected cluster(s)
    if isempty(obj.selected) || ~obj.hasFig('FigWav')
        return;
    end

    hFigWav = obj.hFigs('FigWav');
    if numel(obj.selected) == 1 % we've selected just one cluster (primary)
        hFigWav.rmPlot('selected2'); % if already selected, hide it
        plotSelectedMeans(hFigWav, obj.hClust, obj.selected, 'selected1', obj.maxAmp, obj.hCfg);
    else % just plot #2 for now
        plotSelectedMeans(hFigWav, obj.hClust, obj.selected(1), 'selected1', obj.maxAmp, obj.hCfg);
        plotSelectedMeans(hFigWav, obj.hClust, obj.selected(2), 'selected2', obj.maxAmp, obj.hCfg);
    end
end

%% LOCAL FUNCTIONS
function iCluster = plotSelectedMeans(hFigWav, hClust, iCluster, plotKey, maxAmp, hCfg)
    %PLOTSELECTEDMEANS Plot an overlay on selected cluster's mean traces
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

    hFigWav.multiplot(plotKey, maxAmp, jrclust.views.getXRange(iCluster, hCfg), meanWf);
    hFigWav.plotApply(plotKey, @uistack, 'top');
end