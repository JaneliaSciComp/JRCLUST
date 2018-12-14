function hFigISI = doPlotFigISI(hFigISI, hClust, hCfg, selected)
    %DOPLOTFIGISI Plot return map (ISI vs. previous ISI)
    if numel(selected) == 1
        iCluster = selected(1);
        jCluster = iCluster;
    else
        iCluster = selected(1);
        jCluster = selected(2);
    end

    [iIsiK, iIsiK1] = getReturnMap(iCluster, hClust, hCfg);
    if isempty(hFigISI.figData)
        hFigISI.axes();
        %S_fig.hPlot1 = plot(S_fig.hAx, nan, nan, 'ko');
        hFigISI.addPlot('hPlot1', nan, nan, 'ko');
        %S_fig.hPlot2 = plot(S_fig.hAx, nan, nan, 'ro');
        hFigISI.addPlot('hPlot2', nan, nan, 'ro');

        %set(S_fig.hAx, 'XScale','log', 'YScale','log');
        hFigISI.axSet('XScale','log', 'YScale','log');

        %xlabel('ISI_{k} (ms)'); ylabel('ISI_{k+1} (ms)');
        hFigISI.xlabel('ISI_{k} (ms)');
        hFigISI.ylabel('ISI_{k+1} (ms)');

        %axis_(S_fig.hAx, [1 10000 1 10000]);
        hFigISI.axis([1 10000 1 10000]);
        hFigISI.grid('on');

        % show refractory line
        %line(get(S_fig.hAx,'XLim'), hCfg.spkRefrac_ms*[1 1], 'Color', [1 0 0]);
        hFigISI.addLine('refracLine1', hFigISI.axGet('XLim'), hCfg.refracIntms*[1 1], 'Color', [1 0 0]);
        %line(hCfg.spkRefrac_ms*[1 1], get(S_fig.hAx,'YLim'), 'Color', [1 0 0]);
        hFigISI.addLine('refracLine2', hCfg.refracIntms*[1 1], hFigISI.axGet('YLim'), 'Color', [1 0 0]);
    end

    hFigISI.updatePlot('hPlot1', iIsiK, iIsiK1);

    if iCluster ~= jCluster
        [jIsiK, jIsiK1] = getReturnMap(jCluster, hClust, hCfg);
        hFigISI.updatePlot('hPlot2', jIsiK, jIsiK1);
    else
        hFigISI.hidePlot('hPlot2');
    end

    %set(hFigISI, 'UserData', S_fig);
end

%% LOCAL FUNCTIONS
function [isiK, isiK1] = getReturnMap(iCluster, hClust, hCfg)
    %GETRETURNMAP subset 
    clusterTimes = double(hClust.spikeTimes(hClust.spikesByCluster{iCluster}))/hCfg.sampleRate;
    clusterISIMs = diff(clusterTimes * 1000); % in msec

    %isiK = clusterISIMs(1:end-1);
    %isiK1 = clusterISIMs(2:end);
    %subset = randperm(numel(isiK), min(hCfg.nShow, numel(isiK)));
    subset = randperm(numel(clusterISIMs) - 1, min(hCfg.nShow, numel(clusterISIMs) - 1));
    %isiK = isiK(subset);
    isiK = clusterISIMs(subset);
    %isiK1 = isiK1(subset);
    isiK1 = clusterISIMs(subset+1);
end
