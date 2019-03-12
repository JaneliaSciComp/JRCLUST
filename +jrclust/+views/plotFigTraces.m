function tracesFilt = plotFigTraces(hFigTraces, hCfg, tracesRaw, resetAxis, hClust)
    %PLOTFIGTRACES Plot raw traces view
    hBox = jrclust.utils.qMsgBox('Plotting...', 0, 1);

    hFigTraces.wait(1);

    sampleRate = hCfg.sampleRate / hCfg.nSkip;
    viSamples1 = 1:hCfg.nSkip:size(tracesRaw, 2);
    evtWindowSamp = round(hCfg.evtWindowSamp / hCfg.nSkip); %show 2x of range

    if strcmpi(hFigTraces.figData.filter, 'on')
        % back up old settings
        sampleRateOld = hCfg.sampleRate;
        filterTypeOld = hCfg.filterType;

        % temporarily alter settings to send traces through a filter
        hCfg.sampleRate = sampleRate;
        hCfg.useGPU = 0;
        hCfg.filterType = hCfg.dispFilter;

        if hCfg.fftThresh > 0
            tracesRaw = jrclust.filters.fftClean(tracesRaw, hCft.fftThresh, hCfg);
        end

        tracesFilt = jrclust.filters.filtCAR(tracesRaw(:, viSamples1), [], [], 0, hCfg);
        tracesFilt = jrclust.utils.bit2uV(tracesFilt, hCfg);
        filterToggle = hCfg.filterType;

        % restore old settings
        hCfg.sampleRate = sampleRateOld;
        hCfg.useGPU = 1;
        hCfg.filterType = filterTypeOld;
    else
        tracesFilt = jrclust.utils.meanSubtract(single(tracesRaw(:, viSamples1))) * hCfg.bitScaling;
        filterToggle = 'off';
    end

    if hCfg.nSegmentsTraces == 1
        XData = ((hFigTraces.figData.windowBounds(1):hCfg.nSkip:hFigTraces.figData.windowBounds(end))-1) / hCfg.sampleRate;
        XLabel = 'Time (s)';
    else
        XData = (0:(size(tracesFilt, 2) - 1)) / (hCfg.sampleRate / hCfg.nSkip) + (hFigTraces.figData.windowBounds(1)-1) / hCfg.sampleRate;
        [multiBounds, multiRange, multiEdges] = jrclust.views.sampleSkip(hFigTraces.figData.windowBounds, hFigTraces.figData.nSamplesTotal, hCfg.nSegmentsTraces);

        tlim_show = (cellfun(@(x) x(1), multiBounds([1, end]))) / hCfg.sampleRate;
        XLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', hCfg.nSegmentsTraces, tlim_show, diff(hCfg.dispTimeLimits));

        mrX_edges = XData(repmat(multiEdges(:)', [3, 1]));
        mrY_edges = repmat([0; hCfg.nSites + 1; nan], 1, numel(multiEdges));

        hFigTraces.plotApply('hEdges', @set, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
        csTime_bin = cellfun(@(x) sprintf('%0.1f', x(1)/hCfg.sampleRate), multiBounds, 'UniformOutput', 0);
        hFigTraces.axApply('default', @set, {'XTick', 'XTickLabel'}, {XData(multiEdges), csTime_bin});
    end

    hFigTraces.multiplot('hPlot', hFigTraces.figData.maxAmp, XData, tracesFilt', 1:hCfg.nSites);

    hFigTraces.axApply('default', @grid, hFigTraces.figData.grid);
    hFigTraces.axApply('default', @set, 'YTick', 1:hCfg.nSites);
    hFigTraces.axApply('default', @title, sprintf(hFigTraces.figData.title, hFigTraces.figData.maxAmp));
    hFigTraces.axApply('default', @xlabel, XLabel);
    hFigTraces.axApply('default', @ylabel, 'Site #');
    hFigTraces.plotApply('hPlot', @set, 'Visible', hFigTraces.figData.traces);

    % Delete spikes from other threads (TODO: break this out into a function)
    plotKeys = keys(hFigTraces.hPlots);
    chSpk = plotKeys(startsWith(plotKeys, 'chSpk'));
    if ~isempty(chSpk)
        cellfun(@(pk) hFigTraces.rmPlot(pk), chSpk);
    end

    % plot spikes
    if strcmpi(hFigTraces.figData.spikes, 'on') && ~isempty(hClust)
        recPos = find(strcmp(hFigTraces.figData.hRec.rawPath, hCfg.rawRecordings));
        if recPos == 1
            offset = 0;
        else
            % find all recordings coming before hRec and sum up nSamples
            % for each
            hRecs = arrayfun(@(iRec) jrclust.detect.newRecording(hCfg.rawRecordings{iRec}, hCfg), 1:(recPos-1), 'UniformOutput', 0);
            offset = sum(cellfun(@(hR) hR.nSamples, hRecs));
        end

        recTimes = hClust.spikeTimes - int32(offset);

        tStart = single(hFigTraces.figData.windowBounds(1) - 1)/hCfg.sampleRate;
        if hCfg.nSegmentsTraces > 1
            spikesInRange = inRange(recTimes, multiBounds);
            spikeSites = hClust.spikeSites(spikesInRange);
            spikeTimes = double(recTimes(spikesInRange));
            spikeTimes = round(whereMember(spikeTimes, multiRange) / hCfg.nSkip);
        else
            spikesInRange = recTimes >= hFigTraces.figData.windowBounds(1) & recTimes < hFigTraces.figData.windowBounds(end);
            spikeSites = hClust.spikeSites(spikesInRange);
            spikeTimes = double(recTimes(spikesInRange));
            spikeTimes = round((spikeTimes - hCfg.sampleRate*tStart) / hCfg.nSkip); % time offset
        end

        spikeSites = single(spikeSites);

        % check if clustered
        if isempty(hClust)
            for iSite = 1:hCfg.nSites % deal with subsample factor
                onSite = find(spikeSites == iSite);
                if isempty(onSite)
                    continue;
                end

                timesOnSite = spikeTimes(onSite);
                [mrY11, mrX11] = vr2mr3_(tracesFilt(iSite, :), timesOnSite, evtWindowSamp); %display purpose x2

                mrT11 = single(mrX11-1) / sampleRate + tStart;

                plotKey = sprintf('chSpk%d', iSite);
                hFigTraces.addPlot(plotKey, @line, ...
                                   nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5);

                hFigTraces.multiplot(plotKey, hFigTraces.figData.maxAmp, mrT11, mrY11, iSite);
            end
        else % different color for each cluster
            inRangeClusters = hClust.spikeClusters(spikesInRange);
            spikeColors = [jet(hClust.nClusters); 0 0 0];
            lineWidths = (mod((1:hClust.nClusters) - 1, 3) + 1)'/2 + 0.5;  %(randi(3, S_clu.nClusters, 1)+1)/2;

            % shuffle colors
            spikeColors = spikeColors(randperm(size(spikeColors, 1)), :);
            lineWidths = lineWidths(randperm(size(lineWidths, 1)), :);

            nSpikes = numel(spikeTimes);

            for iSpike = 1:nSpikes
                iCluster = inRangeClusters(iSpike);
                if iCluster <= 0
                    continue;
                end

                iTime = spikeTimes(iSpike);
                iSite = spikeSites(iSpike);
                iColor = spikeColors(iCluster, :);
                iLinewidth = lineWidths(iCluster);

                [mrY11, mrX11] = vr2mr3_(tracesFilt(iSite, :), iTime, evtWindowSamp); %display purpose x2
                mrT11 = double(mrX11-1) / sampleRate + tStart;

                plotKey = sprintf('chSpk%d', iSpike);
                hFigTraces.addPlot(plotKey, @line, ...
                                   nan, nan, 'Color', iColor, 'LineWidth', iLinewidth);
                hFigTraces.multiplot(plotKey, hFigTraces.figData.maxAmp, mrT11, mrY11, iSite);
            end
        end
    end

    if resetAxis
        jrclust.views.resetFigTraces(hFigTraces, tracesRaw, hCfg);
    end

    hFigTraces.figApply(@set, 'Name', sprintf('%s: filter: %s', hCfg.configFile, filterToggle));
    hFigTraces.wait(0);

    jrclust.utils.tryClose(hBox);
end

%% LOCAL FUNCTIONS
function isInRange = inRange(vals, ranges)
    isInRange = false(size(vals));

    if ~iscell(ranges)
        ranges = {ranges};
    end

    for iRange = 1:numel(ranges)
        bounds = ranges{iRange};
        isInRange = isInRange | (vals >= bounds(1) & vals <= bounds(2));
    end
end

function loc = whereMember(needle, haystack)
    % needle = int32(needle); % ??
    % haystack = int32(haystack); % ??
    [~, loc] = ismember(int32(needle), int32(haystack));
end

function [mr, ranges] = vr2mr3_(traces, spikeTimes, evtWindow)
    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype

    % prepare indices
    spikeTimes = spikeTimes(:)';
    ranges = int32(bsxfun(@plus, (evtWindow(1):evtWindow(end))', spikeTimes));

    ranges(ranges < 1) = 1;
    ranges(ranges > numel(traces)) = numel(traces); %keep # sites consistent

    % build spike table
    nSpikes = numel(spikeTimes);
    mr = reshape(traces(ranges(:)), [], nSpikes);
end
