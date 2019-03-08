function updateFigISI(obj)
    %UPDATEFIGISI Plot return map
    if isempty(obj.selected) || ~obj.hasFig('FigISI')
        return;
    end

    plotFigISI(obj.hFigs('FigISI'), obj.hClust, obj.hCfg, obj.selected);
end

%% LOCAL FUNCTIONS
function hFigISI = plotFigISI(hFigISI, hClust, hCfg, selected)
    %DOPLOTFIGISI Plot return map (ISI vs. previous ISI)
    if numel(selected) == 1
        iCluster = selected(1);
        jCluster = iCluster;
    else
        iCluster = selected(1);
        jCluster = selected(2);
    end

    [iIsiK, iIsiK1] = getReturnMap(iCluster, hClust, hCfg);
    if ~hFigISI.hasAxes('default')
        hFigISI.addAxes('default');
        hFigISI.addPlot('foreground', nan, nan, 'Color', hCfg.colorMap(2, :), 'Marker', 'o', 'LineStyle', 'none');
        hFigISI.addPlot('foreground2', nan, nan, 'Color', hCfg.colorMap(3, :), 'Marker', 'o', 'LineStyle', 'none');

        hFigISI.axApply('default', @set, 'XScale','log', 'YScale','log');

        hFigISI.axApply('default', @xlabel, 'ISI_{k} (ms)');
        hFigISI.axApply('default', @ylabel, 'ISI_{k+1} (ms)');

        %axis_(S_fig.hAx, [1 10000 1 10000]);
        hFigISI.axApply('default', @axis, [1 10000 1 10000]);
        hFigISI.axApply('default', @grid, 'on');

        % show refractory line
        hFigISI.addPlot('refracLine1', @line, hFigISI.axApply('default', @get, 'XLim'), hCfg.refracInt*[1 1], 'Color', [1 0 0]);
        hFigISI.addPlot('refracLine2', @line, hCfg.refracInt*[1 1], hFigISI.axApply('default', @get, 'YLim'), 'Color', [1 0 0]);
    end

    hFigISI.updatePlot('foreground', iIsiK, iIsiK1);

    if iCluster ~= jCluster
        [jIsiK, jIsiK1] = getReturnMap(jCluster, hClust, hCfg);
        hFigISI.updatePlot('foreground2', jIsiK, jIsiK1);
    else
        hFigISI.clearPlot('foreground2');
    end
end

function [isiK, isiK1] = getReturnMap(iCluster, hClust, hCfg)
    %GETRETURNMAP subset 
    clusterTimes = double(hClust.spikeTimes(hClust.spikesByCluster{iCluster}))/hCfg.sampleRate;
    clusterISI = diff(clusterTimes * 1000); % in msec

    %isiK = clusterISIMs(1:end-1);
    %isiK1 = clusterISIMs(2:end);
    %subset = randperm(numel(isiK), min(hCfg.nShow, numel(isiK)));
    subset = randperm(numel(clusterISI) - 1, min(hCfg.nSpikesFigISI, numel(clusterISI) - 1));
    %isiK = isiK(subset);
    isiK = clusterISI(subset);
    %isiK1 = isiK1(subset);
    isiK1 = clusterISI(subset+1);
end
