function autoSplit(obj, multisite)
    %AUTOSPLIT
    if numel(obj.selected) > 1
        return;
    end
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    end

    if obj.hClust.unitCount(obj.selected) < 2
        msgbox('At least two spikes required for splitting');
        return;
    end

    iCluster = obj.selected(1);
    iSite = obj.hClust.clusterSites(iCluster);

    if multisite
        spikeSites = obj.hCfg.siteNeighbors(1:end - obj.hCfg.nSitesExcl, iSite);
    else
        spikeSites = iSite;
    end

    iSpikes = obj.hClust.spikesByCluster{iCluster};

    sampledTraces = obj.hClust.getSpikeWindows(iSpikes, spikeSites, 0, 1);
    sampledTraces = reshape(sampledTraces, [], size(sampledTraces, 3));

    % get Vpp of cluster spikes on current site (in FigTime)
    localSpikes = squeeze(obj.hClust.getSpikeWindows(iSpikes, obj.currentSite, 0, 1));
    localVpp = max(localSpikes) - min(localSpikes); % TW calculate amplitudes on the fly

    clusterTimes = obj.hClust.spikeTimes(iSpikes);

    % create split figure
    hFigSplit = jrclust.views.Figure('FigSplit', [0.05 0.05 0.8 0.8], sprintf('Split cluster %d', iCluster), 0, 0);
    hFigSplit.figApply(@set, 'Visible', 'off');
    hFigSplit.figData.nSplits = 0;

    % create a dialog
    splitDlg = dialog('Name', 'Split cluster', ...
                      'Units', 'Normalized', ...
                      'Position', [0.4, 0.5, 0.2, 0.15]);

    uicontrol('Parent', splitDlg, 'Style', 'text', ...
              'String', 'Number of splits', ...
              'Units', 'Normalized', ...
              'Position', [0.4 0.4 0.2 0.07]);

    nsplit = uicontrol('Parent', splitDlg, 'Style', 'edit', ...
                       'String', 2, ...
                       'Units', 'Normalized', ...
                       'Position', [0.4 0.25 0.2 0.15]);

    % button group: K-means, K-medoids, Hierarchical clustering
    btngrp = uibuttongroup('Parent', splitDlg, 'Units', 'Normalized', ...
                           'Position', [0.1 0.1 0.3 0.8]);
    uicontrol('Parent', btngrp, 'Style', 'radiobutton', 'String', 'Hierarchical (Ward)', ...
              'Units', 'Normalized', 'Position', [0.1, 0.1, 1, 0.15]);
    uicontrol('Parent', btngrp, 'Style', 'radiobutton', 'String', 'K-means', ...
              'Units', 'Normalized', 'Position', [0.1, 0.7, 1, 0.15]);
    uicontrol('Parent', btngrp, 'Style', 'radiobutton', 'String', 'K-medoids', ...
              'Units', 'Normalized', 'Position', [0.1, 0.4, 1, 0.15]);

    uicontrol('Parent', splitDlg, 'Style', 'pushbutton', ...
              'String', 'Split', ...
              'Units', 'Normalized', ...
              'Position', [0.4 0.1 0.2 0.15], ...
              'Callback', @(hO, hE) preSplit(splitDlg, nsplit.String, btngrp, hFigSplit));
    uicontrol('Parent', splitDlg, 'Style', 'pushbutton', ...
              'String', 'Cancel', ...
              'Units', 'Normalized', ...
              'Position', [0.6 0.1 0.2 0.15], ...
              'Callback', @(hO, hE) doCancel(splitDlg, hFigSplit));
    uicontrol('Parent', splitDlg, 'Style', 'text', ...
              'String', 'How to recluster this unit?', ...
              'Units', 'Normalized', ...
              'Position', [0.1 0.9 0.3 0.07]);

    uiwait(splitDlg);

    if hFigSplit.figData.nSplits <= 1
        return;
    else
        nSplits = hFigSplit.figData.nSplits;
    end

    hBox = jrclust.utils.qMsgBox('Splitting... (this closes automatically)', 0, 1);

    % add plots
    hFigSplit.addAxes('vppTime', 'Units', 'Normalized', 'OuterPosition', [0 0 0.7 0.5]);
    hFigSplit.addAxes('pc12', 'Units', 'Normalized', 'OuterPosition', [0 0.5 0.2333 0.5]);
    hFigSplit.addAxes('pc13', 'Units', 'Normalized', 'OuterPosition', [0.2333 0.5 0.2333 0.5]);
    hFigSplit.addAxes('pc23', 'Units', 'Normalized', 'OuterPosition', [0.4667 0.5 0.2333 0.5]);
    hFigSplit.addAxes('isi', 'Units', 'Normalized', 'OuterPosition', [0.7 0.65 0.3 0.35]);
    hFigSplit.addAxes('meanwf', 'Units', 'Normalized', 'OuterPosition', [0.7 0.4 0.3 0.25]);

    % add controls
    hFigSplit.figData.clustList = hFigSplit.figApply(@uicontrol, 'Style', 'listbox', ...
                                                     'Units', 'normalized', ...
                                                     'Position', [0.7 0.05 0.1 0.25], ...
                                                     'String', num2str((1:nSplits)'), ...
                                                     'Value', 1, 'BackgroundColor', [1 1 1], ...
                                                     'Max', 3, 'Min', 1, ...
                                                     'Callback',  @(hO, hE) updateSplitPlots(hFigSplit));

                                                 
    hFigSplit.figData.doneButton = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                      'Units', 'normalized', ...
                                                      'Position', [0.84 0.25 0.12 0.05], ...
                                                      'String', 'Keep splits', ...
                                                      'FontWeight','Normal', ...
                                                      'Callback', @(hO, hE) finishSplit(hFigSplit, 0));
    hFigSplit.figData.quitButton = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                      'Units', 'normalized', ...
                                                      'Position', [0.84 0.2 0.12 0.05], ...
                                                      'String', 'Discard splits', ...
                                                      'FontWeight','Normal', ...
                                                      'Callback', @(hO, hE) finishSplit(hFigSplit, 1));

    hFigSplit.figData.mergeButton = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                       'Units', 'normalized', ...
                                                       'Position', [0.84 0.1 0.12 0.05], ...
                                                       'String', 'Merge selected', ...
                                                       'FontWeight','Normal', ...
                                                       'Callback', @(hO, hE) mergeSelected(hFigSplit));
    hFigSplit.figData.button12 = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                    'Units', 'normalized', ...
                                                    'Position', [0.84 0.05 0.04 0.05], ...
                                                    'String', 'Manual 2 vs. 1', ...
                                                    'FontWeight','Normal', ...
                                                    'Callback', @(hO, hE) manualSplit(hFigSplit, [1 2]));
    hFigSplit.figData.button13 = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                    'Units', 'normalized', ...
                                                    'Position', [0.88 0.05 0.04 0.05], ...
                                                    'String', 'Manual 3 vs. 1', ...
                                                    'FontWeight','Normal', ...
                                                    'Callback', @(hO, hE) manualSplit(hFigSplit, [1 3]));
    hFigSplit.figData.button23 = hFigSplit.figApply(@uicontrol, 'Style', 'pushbutton', ...
                                                    'Units', 'normalized', ...
                                                    'Position', [0.92 0.05 0.04 0.05], ...
                                                    'String', 'Manual 3 vs. 2', ...
                                                    'FontWeight','Normal', ...
                                                    'Callback', @(hO, hE) manualSplit(hFigSplit, [2 3]));

    hFigSplit.figData.finished = 0;
    hFigSplit.figData.iSite = iSite;
    hFigSplit.figData.refracInt = obj.hCfg.refracInt/1000;
    while 1
        [assigns, pcaFeatures] = doAutoSplit(sampledTraces, [double(clusterTimes) double(localVpp')], hFigSplit); %TW
        hFigSplit.figData.pcaFeatures = pcaFeatures;
        hFigSplit.figData.clusterTimes = double(clusterTimes)/obj.hCfg.sampleRate;
        hFigSplit.figData.localVpp = double(localVpp');
        hFigSplit.figData.assignPart = arrayfun(@(i) find(assigns == i)', 1:max(assigns), 'UniformOutput', 0);
        hFigSplit.figData.sampledTraces = localSpikes;

        jrclust.utils.tryClose(hBox);

        updateSplitPlots(hFigSplit);
        hFigSplit.figApply(@uiwait);
        if ~hFigSplit.isReady % user closed the window
            assignPart = [];
            break;
        elseif hFigSplit.figData.finished
            assignPart = hFigSplit.figData.assignPart; 
            break;
        end
    end

    hFigSplit.close();

    obj.isWorking = 0;
    if ~isempty(assignPart)
        obj.splitCluster(iCluster, assignPart);
    end
end

%% LOCAL FUNCTIONS
function [assigns, pcaFeatures] = doAutoSplit(sampledSpikes, spikeFeatures, hFigSplit)
    %DOAUTOSPLIT
    % TODO: ask users number of clusters and split multi-way
    %Make automatic split of clusters using PCA + hierarchical clustering
    [~, pcaFeatures, ~] = pca(double(sampledSpikes'), 'Centered', 1, 'NumComponents', 3);
    combinedFeatures = [spikeFeatures pcaFeatures];
    nSpikes = size(sampledSpikes, 2);

    % ask how many clusters there are
    try
        combinedFeatures = (combinedFeatures - mean(combinedFeatures, 1)) ./ std(combinedFeatures, 1);
        assigns = hFigSplit.figData.hFunSplit(combinedFeatures);
    catch % not enough features to automatically split or some other failure
        assigns = ones(nSpikes, 1);
    end
end

function doCancel(splitDlg, hFigSplit)
    delete(splitDlg);
    hFigSplit.figData.nSplits = 0;
end

function preSplit(splitDlg, nsplit, btngrp, hFigSplit)
    whichSelected = logical(arrayfun(@(child) child.Value, btngrp.Children));
    selected = btngrp.Children(whichSelected);

    nsplit = floor(str2double(nsplit));
    if isnan(nsplit)
        nsplit = 2;
    end

    switch selected.String
        case 'K-means'
            hFigSplit.figData.hFunSplit = @(X) kmeans(X, nsplit);
        case 'K-medoids'
            hFigSplit.figData.hFunSplit = @(X) kmedoids(X, nsplit);
        otherwise
            hFigSplit.figData.hFunSplit = @(X) clusterdata(X, ...
                                                           'maxclust', nsplit, ...
                                                           'linkage', 'ward', ...
                                                           'distance', 'euclidean', ...
                                                           'savememory', 'on');
    end

    delete(splitDlg)

    hFigSplit.figData.nSplits = nsplit;
end

function updateSplitPlots(hFigSplit)
    clustList = hFigSplit.figData.clustList;

    assignPart = hFigSplit.figData.assignPart;
    clusterTimes = hFigSplit.figData.clusterTimes;
    localVpp = hFigSplit.figData.localVpp;
    pcaFeatures = hFigSplit.figData.pcaFeatures;
    sampledTraces = hFigSplit.figData.sampledTraces;
    refracInt = hFigSplit.figData.refracInt;

    nSplits = numel(assignPart);

    selectedClusters = get(clustList, 'Value');
    selectedSpikes = sort([assignPart{selectedClusters}]);

    cmap = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.3010, 0.7450, 0.9330; 0.6350, 0.0780, 0.1840];
    legends = cell(nSplits, 1);

    % plot split features
    hFigSplit.axApply('vppTime', @cla);
    for i = 1:nSplits
        partSpikes = assignPart{i};
        skip = ceil(numel(partSpikes) / 1e4);

        XData = clusterTimes(partSpikes);
        YData = localVpp(partSpikes);
        hFigSplit.addPlot(sprintf('plot1%d', i), @plot, ...
                          hFigSplit.hAxes('vppTime'), ...
                          XData(1:skip:end), YData(1:skip:end), ...
                          '.', 'Color', cmap(mod(i-1, size(cmap, 1)) + 1, :), ...
                          'MarkerSize', jrclust.utils.ifEq(ismember(i, selectedClusters), 10, 5));

        hFigSplit.axApply('vppTime', @hold, 'on');

        legends{i} = sprintf('Split %d (%d spikes)', i, numel(partSpikes));
    end

    hFigSplit.axApply('vppTime', @legend, legends, 'Location', 'Best');
    hFigSplit.axApply('vppTime', @set, 'YDir', 'Reverse');
    hFigSplit.axApply('vppTime', @view, 0, 90);
    hFigSplit.axApply('vppTime', @xlim, [min(clusterTimes) max(clusterTimes)]);
    hFigSplit.axApply('vppTime', @ylim, [min(localVpp) max(localVpp)]);
    hFigSplit.axApply('vppTime', @xlabel, 'Spike times (s)');
    hFigSplit.axApply('vppTime', @ylabel, 'Vpp (\muV)');

    hFigSplit.axApply('pc12', @cla);
    hFigSplit.axApply('pc13', @cla);
    hFigSplit.axApply('pc23', @cla);
    for i = 1:nSplits
        partSpikes = assignPart{i};
        skip = ceil(numel(partSpikes) / 1e4);

        % PC1vs2
        XData = pcaFeatures(partSpikes, 2);
        YData = pcaFeatures(partSpikes, 1);
        hFigSplit.addPlot(sprintf('pc12-%d', i), @plot, ...
                          hFigSplit.hAxes('pc12'), ...
                          XData(1:skip:end), YData(1:skip:end), ...
                          '.', 'Color', cmap(mod(i-1, size(cmap, 1)) + 1, :), ...
                          'MarkerSize', jrclust.utils.ifEq(ismember(i, selectedClusters), 10, 5));

        hFigSplit.axApply('pc12', @hold, 'on');

        % PC1vs3
        XData = pcaFeatures(partSpikes, 3);
        YData = pcaFeatures(partSpikes, 1);
        hFigSplit.addPlot(sprintf('pc13-%d', i), @plot, ...
                          hFigSplit.hAxes('pc13'), ...
                          XData(1:skip:end), YData(1:skip:end), ...
                          '.', 'Color', cmap(mod(i-1, size(cmap, 1)) + 1, :), ...
                          'MarkerSize', jrclust.utils.ifEq(ismember(i, selectedClusters), 10, 5));

        hFigSplit.axApply('pc13', @hold, 'on');

        % PC2vs3
        XData = pcaFeatures(partSpikes, 3);
        YData = pcaFeatures(partSpikes, 2);
        hFigSplit.addPlot(sprintf('pc23-%d', i), @plot, ...
                          hFigSplit.hAxes('pc23'), ...
                          XData(1:skip:end), YData(1:skip:end), ...
                          '.', 'Color', cmap(mod(i-1, size(cmap, 1)) + 1, :), ...
                          'MarkerSize', jrclust.utils.ifEq(ismember(i, selectedClusters), 10, 5));

        hFigSplit.axApply('pc23', @hold, 'on');
    end

    hFigSplit.axApply('pc12', @hold, 'off');
    hFigSplit.axApply('pc12', @xlabel, 'PC 2 (AU)');
    hFigSplit.axApply('pc12', @ylabel, 'PC 1 (AU)');

    hFigSplit.axApply('pc13', @hold, 'off');
    hFigSplit.axApply('pc13', @xlabel, 'PC 3 (AU)');
    hFigSplit.axApply('pc13', @ylabel, 'PC 1 (AU)');

    hFigSplit.axApply('pc23', @hold, 'off');
    hFigSplit.axApply('pc23', @xlabel, 'PC 3 (AU)');
    hFigSplit.axApply('pc23', @ylabel, 'PC 2 (AU)');

    % plot ISI
    selectedTimes = clusterTimes(selectedSpikes);
    selectedISI = diff([selectedTimes; inf]);

    if isempty(selectedISI)
        selectedISI = 1;
    end

    ISIedges = -0.02:0.0005:0.02;
    nISI = histc([selectedISI; -selectedISI], ISIedges);

    mx = 1.05 .* max(nISI);

    if ~hFigSplit.hasPlot('isiHist')
        % plot ISI region (region in hCfg.refracInt)
        hFigSplit.addPlot('isiRegion', @fill, ...
                      hFigSplit.hAxes('isi'), ...
                      [-refracInt refracInt refracInt -refracInt], [0 0 mx mx], 'r');
        hFigSplit.plotApply('isiRegion', @set, 'LineStyle', 'none');

        % plot histogram of ISIs in cluster
        hFigSplit.addPlot('isiHist', @bar, ...
                          hFigSplit.hAxes('isi'), ISIedges + mean(diff(ISIedges))./2, ...
                          nISI, 'BarWidth', 1, 'FaceAlpha', 0.75);
        hFigSplit.axApply('isi', @xlim, [-0.02 0.02]);
        hFigSplit.axApply('isi', @xlabel, 'ISI (s)');
        hFigSplit.axApply('isi', @ylabel, 'Counts');
        hFigSplit.axApply('isi', @title, 'ISI histogram (refractory interval in red)');
    else
        hFigSplit.updatePlot('isiRegion', [-refracInt refracInt refracInt -refracInt], [0 0 mx mx]);
        hFigSplit.updatePlot('isiHist', ISIedges + mean(diff(ISIedges))./2, nISI);
    end

    violations = false(size(selectedISI));
    violations(1:end-1) = selectedISI(1:end-1) < refracInt;
    violations(2:end)   = violations(2:end) | selectedISI(1:end - 1) < refracInt;

    violAll = false(size(selectedSpikes));
    violAll(selectedSpikes) = violations;

    % highlight violations on feature plot
    if any(violAll)
        violAll = find(violAll);
        hFigSplit.addPlot(sprintf('plot1%d-viol', i), @plot, ...
                          hFigSplit.hAxes('vppTime'), ...
                          clusterTimes(violAll), ...
                          localVpp(violAll), ...
                          'ko', 'LineWidth', 2);
        legends{end+1} = 'ISI violations';
    end

    hFigSplit.axApply('vppTime', @hold, 'off');
    hFigSplit.axApply('vppTime', @legend, legends, 'Location', 'Best');

    % plot mean traces
    unselectedClusters = setdiff(1:nSplits, selectedClusters);
    hFigSplit.axApply('meanwf', @cla);
    nSpikesPlot = 50;

    % plot unselected traces first...
    for ii = 1:numel(unselectedClusters)
        i = unselectedClusters(ii);
        partSpikes = assignPart{i};
        hFigSplit.addPlot(sprintf('traces-%d', i), @plot, ...
                          hFigSplit.hAxes('meanwf'), ...
                          sampledTraces(:, randsample(partSpikes, min(numel(partSpikes), nSpikesPlot))), ...
                          'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);

        hFigSplit.axApply('meanwf', @hold, 'on');
    end

    % ...plot selected traces on top of these...
    for ii = 1:numel(selectedClusters)
        i = selectedClusters(ii);
        iColor = cmap(mod(i-1, size(cmap,1))+1, :)+(1-cmap(mod(i-1, size(cmap,1))+1, :))./2;

        partSpikes = assignPart{i};
        hFigSplit.addPlot(sprintf('traces-%d', i), @plot, ...
                          hFigSplit.hAxes('meanwf'), ...
                          sampledTraces(:, randsample(partSpikes, min(numel(partSpikes), nSpikesPlot))), ...
                          'Color', iColor, 'LineWidth', 0.5);

        hFigSplit.axApply('meanwf', @hold, 'on');
    end

    % ...plot unselected mean waveforms here...
    for ii = 1:numel(unselectedClusters)
        i = unselectedClusters(ii);
        partSpikes = assignPart{i};

        hFigSplit.addPlot(sprintf('mean-%d', i), @plot, ...
                          hFigSplit.hAxes('meanwf'), ...
                          mean(sampledTraces(:, partSpikes), 2), ...
                          'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    end

    for ii = 1:numel(selectedClusters)
        i = selectedClusters(ii);
        iColor = cmap(mod(i-1, size(cmap, 1))+1, :);
        partSpikes = assignPart{i};

        hFigSplit.addPlot(sprintf('mean-%d', i), @plot, ...
                          hFigSplit.hAxes('meanwf'), ...
                          mean(sampledTraces(:, partSpikes), 2), ...
                          'Color', iColor, 'LineWidth', 2);
    end

    hFigSplit.axApply('meanwf', @hold, 'off');
    hFigSplit.axApply('meanwf', @xlabel, 'Sample');
    hFigSplit.axApply('meanwf', @ylabel, 'Amplitude (\muV)');
    hFigSplit.axApply('meanwf', @title, sprintf('Mean waveform on site %d', hFigSplit.figData.iSite));
end

function manualSplit(hFigSplit, pcPair)
    clusterTimes = hFigSplit.figData.clusterTimes;
    localVpp = hFigSplit.figData.localVpp;
    pcaFeatures = hFigSplit.figData.pcaFeatures;

    % clear out all axes
    hFigSplit.axApply('pc12', @cla);
    hFigSplit.axApply('pc13', @cla);
    hFigSplit.axApply('pc23', @cla);
    hFigSplit.axApply('vppTime', @cla);
    hFigSplit.axApply('isi', @cla);
    hFigSplit.axApply('meanwf', @cla);

    nSpikes = numel(clusterTimes);
    skip = ceil(numel(nSpikes) / 1e4);

    % PC1vs2
    XData = pcaFeatures(:, 2);
    YData = pcaFeatures(:, 1);
    hFigSplit.addPlot('pc12', @plot, ...
                      hFigSplit.hAxes('pc12'), ...
                      XData(1:skip:end), YData(1:skip:end), ...
                      '.', 'Color', 'k');

    hFigSplit.axApply('pc12', @hold, 'off');

    % PC1vs3
    XData = pcaFeatures(:, 3);
    YData = pcaFeatures(:, 1);
    hFigSplit.addPlot('pc13', @plot, ...
                      hFigSplit.hAxes('pc13'), ...
                      XData(1:skip:end), YData(1:skip:end), ...
                      '.', 'Color', 'k');

    hFigSplit.axApply('pc13', @hold, 'off');

    % PC2vs3
    XData = pcaFeatures(:, 3);
    YData = pcaFeatures(:, 2);
    hFigSplit.addPlot('pc23', @plot, ...
                      hFigSplit.hAxes('pc23'), ...
                      XData(1:skip:end), YData(1:skip:end), ...
                      '.', 'Color', 'k');

    hFigSplit.axApply('pc13', @hold, 'off');

    % vpp vs time
    hFigSplit.addPlot('vppTime', @plot, ...
                      hFigSplit.hAxes('vppTime'), ...
                      clusterTimes(1:skip:end), localVpp(1:skip:end), ...
                      '.', 'Color', 'k');

    hFigSplit.axApply('vppTime', @hold, 'off');

    if all(pcPair == [1 2])
        axKey = 'pc12';
        featuresX = pcaFeatures(1:skip:end, 2);
        featuresY = pcaFeatures(1:skip:end, 1);
    elseif all(pcPair == [1 3])
        axKey = 'pc13';
        featuresX = pcaFeatures(1:skip:end, 3);
        featuresY = pcaFeatures(1:skip:end, 1); 
    else % all(pcPair == [2 3])
        axKey = 'pc23';
        featuresX = pcaFeatures(1:skip:end, 3);
        featuresY = pcaFeatures(1:skip:end, 2);
    end

    hFigSplit.axApply(axKey, @hold, 'on');

    % clear axes
    hPoly = hFigSplit.axApply(axKey, @impoly);
    polyPos = getPosition(hPoly);
    try
        delete(hPoly);
    catch
    end

    % partition spikes by polygon
    assigns = inpolygon(featuresX, featuresY, polyPos(:, 1), polyPos(:, 2));
    assignPart = {find(assigns), find(~assigns)};
    assignPartOld = hFigSplit.figData.assignPart;

    % update plots
    nSplits = 2;
    hFigSplit.figData.assignPart = assignPart;
    set(hFigSplit.figData.clustList, 'String', num2str((1:nSplits)'), 'Value', 1);
    updateSplitPlots(hFigSplit);

    dlgAns = questdlg('Keep?', 'Manual split', 'No');
    switch dlgAns
        case {'No', 'Cancel'}
            nSplits = numel(assignPartOld);
            hFigSplit.figData.assignPart = assignPartOld;
            set(hFigSplit.figData.clustList, 'String', num2str((1:nSplits)'), 'Value', 1);
            updateSplitPlots(hFigSplit);
    end

    hFigSplit.figApply(@uiwait);
end

function mergeSelected(hFigSplit)
    clustList = hFigSplit.figData.clustList;
    selectedClusters = get(clustList, 'Value');

    nSplits = numel(hFigSplit.figData.assignPart);

    if numel(selectedClusters) > 1
        unselectedClusters = setdiff(1:nSplits, selectedClusters);
        unmerged = hFigSplit.figData.assignPart(unselectedClusters);
        merged = {sort([hFigSplit.figData.assignPart{selectedClusters}])};

        hFigSplit.figData.assignPart = cat(2, merged, unmerged);

        nSplits = 1 + numel(unmerged);
        set(hFigSplit.figData.clustList, 'String', num2str((1:nSplits)'), 'Value', 1);
        updateSplitPlots(hFigSplit);
    end
end

function finishSplit(hFigSplit, fAbort)
    if fAbort
        % confirm
        dlgAns = questdlg('Really cancel?', 'Confirm', 'No');
        if strcmp(dlgAns, 'Yes')
            hFigSplit.figData.finished = 1;
            hFigSplit.figData.assignPart = [];
            hFigSplit.figApply(@uiresume);
        end
    else
        % confirm
        dlgAns = questdlg('Finished splitting?', 'Confirm', 'No');
        if strcmp(dlgAns, 'Yes')
            hFigSplit.figData.finished = 1;
            hFigSplit.figApply(@uiresume);
        end
    end
end

% function d12 = madDist(features1, features2)
%     % distance between two clusters
%     if ~ismatrix(features1)
%         features1 = reshape(features1, [], size(features1, 3));
%     end
%     if ~ismatrix(features2)
%         features2 = reshape(features2, [], size(features2, 3));
%     end
% 
%     med1 = median(features1, 2);
%     med2 = median(features2, 2);
% 
%     medDiffs = med1 - med2;
%     norm12 = norm(medDiffs, 2);
%     vrFet12_med1 = medDiffs/norm12;
% 
%     mad1 = median(abs(vrFet12_med1' * bsxfun(@minus, features1, med1)));
%     mad2 = median(abs(vrFet12_med1' * bsxfun(@minus, features2, med2)));
% 
%     d12 = norm12 / norm([mad1, mad2], 2);
% end