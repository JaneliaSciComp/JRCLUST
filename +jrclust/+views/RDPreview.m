function hCFig = RDPreview(hClust)
    hCFig = jrclust.views.ControlFigure('FigRLPreview', [0.05 0.0 0.5 1], 'RL Preview', 1, 0);
    set(gcf, 'CloseRequestFcn', @finalize);

    hCFig.addAxes('default');

    % populate the site # dropdown
    hCFig.addUicontrol('siteNoLabel', 'Style', 'text', ...
                       'String', 'Site #', ...
                       'Position', [10, 0, 30, 25]);
    hSite = hCFig.addUicontrol('siteNo', 'Style', 'popup', ...
                               'String', arrayfun(@(n) num2str(n), 1:numel(hClust.spikesBySite), 'UniformOutput', 0), ...
                               'Position', [45, 0, 50, 30], ...
                               'Callback', @setSite);

    % populate the feature projection dropdown
    projections = {'amp vs. time'};
    if strcmp(hClust.hCfg.clusterFeature, 'pca')
        nPCs = hClust.hCfg.nPCsPerSite;
        for i = 1:nPCs
            projections{end+1} = sprintf('PC%d vs. amp', i);
        end
        if nPCs > 1
            for i = 1:nPCs
                for j = i+1:nPCs
                    projections{end+1} = sprintf('PC%d vs. PC%d', i, j);
                end
            end
        end
    end

    hCFig.addUicontrol('projectionLabel', 'Style', 'text', ...
                       'String', 'Projection', ...
                       'Position', [115, 0, 50, 25]);
    hProj = hCFig.addUicontrol('projection', 'Style', 'popup', ...
                               'String', projections, ...
                               'Value', 2, ...
                               'Position', [175, 0, 100, 30], ...
                               'Callback', @setProjection);
    hCFig.axApply('default', @title, 'Projections of spikes onto features');
    hFigRD = jrclust.views.Figure('FigRD', [0.55 0 0.4 0.5], ['Cluster rho-delta: ', hClust.hCfg.sessionName], 0, 0);

    hFigWav = jrclust.views.Figure('FigWav', [0.55 0.5 0.4 0.5], ['Filtered traces: ', hClust.hCfg.sessionName], 0, 0);


    % populate the detrend option dropdown
    detrendOptions = {'none', 'local', 'global', 'regress'};
    hCFig.addUicontrol('detrendLabel', 'Style', 'text', ...
                       'String', 'Detrend mode', ...
                       'Position', [285, 0, 50, 30]);
    currentDetrend = find(strcmp(detrendOptions, hClust.hCfg.RDDetrendMode));
    hDet = hCFig.addUicontrol('projection', 'Style', 'popup', ...
                              'String', detrendOptions, ...
                              'Value', currentDetrend, ...
                              'Position', [345, 0, 100, 30], ...
                              'Callback', @setDetrend);

    % populate the rho slider
    rhoMin = min(log10(hClust.spikeRho));
    rhoMax = max(log10(hClust.spikeRho));
    currentRho = hClust.hCfg.log10RhoCut;

    hCFig.addUicontrol('rhoLabel', 'Style', 'text', ...
                       'String', ['log10 ' char(961) ' cut'], ...
                       'Position', [450, 0, 75, 25]);
    hRho = hCFig.addUicontrol('rhoSlider', 'Style', 'slider', ...
                              'Min', rhoMin, ...
                              'Max', rhoMax, ...
                              'Value', currentRho, ...
                              'Position', [525, 0, 100, 30], ...
                              'Callback', @setRho);

    hRhoVal = hCFig.addUicontrol('rhoVal', 'Style', 'edit', ...
                                 'String', sprintf('%.3f', currentRho), ...
                                 'Position', [625, 0, 50, 30], ...
                                 'Callback', @setRho);

    % populate the delta slider
    deltaMin = min(log10(hClust.spikeDelta));
    deltaMax = max(log10(hClust.spikeDelta));
    currentDelta = hClust.hCfg.log10DeltaCut;

    hCFig.addUicontrol('deltaLabel', 'Style', 'text', ...
                       'String', ['log10 ' char(948) ' cut'], ...
                       'Position', [675, 0, 75, 25]);
    hDelta = hCFig.addUicontrol('deltaSlider', 'Style', 'slider', ...
                                'Min', deltaMin, ...
                                'Max', deltaMax, ...
                                'Value', currentDelta, ...
                                'Position', [750, 0, 100, 30], ...
                                'Callback', @setDelta);
    hDeltaVal = hCFig.addUicontrol('deltaVal', 'Style', 'edit', ...
                                   'String', sprintf('%.3f', currentDelta), ...
                                   'Position', [850, 0, 50, 30], ...
                                   'Callback', @setDelta);

    hSiteClusters = hCFig.addUicontrol('siteClusters', 'Style', 'listbox', ...
                                       'String', '', ...
                                       'Position', [10, 100, 100, 500], ...
                                       'Callback', @setSelectedClusters);

    handles = [hSite, hProj, hDet, hRho, hRhoVal, hDelta, hDeltaVal, hSiteClusters];

    % initial state
    initialSettings = struct('RDDetrendMode', hClust.hCfg.RDDetrendMode, ...
                             'position', size(hClust.history, 1) - 1, ...
                             'rhoCut', hClust.hCfg.log10RhoCut, ...
                             'deltaCut', hClust.hCfg.log10DeltaCut);
    currentProjection = projections{2};
    selectedClusters = [];

    %plotFeatures();
    currentSite = 0;
    ucStr = {};
    setSite(struct('Value', 1));
    plotRD();

%     hCFig.setMouseable();
%     hFigRD.setMouseable();
%     hFigWav.setMouseable();

    function setSite(uic, ~)
        if currentSite == uic.Value
            return;
        end
        currentSite = uic.Value;

        updateSiteClusters();

        plotFeatures();
        plotWaves();
    end

    function setProjection(uic, ~)
        if strcmp(currentProjection, projections{uic.Value})
            return;
        end
        currentProjection = projections{uic.Value};
        plotFeatures();
    end

    function setDetrend(uic, ~)
        if strcmp(currentDetrend, detrendOptions{uic.Value})
            return;
        end
        currentDetrend = detrendOptions{uic.Value};
        hClust.hCfg.RDDetrendMode = currentDetrend;
        doReassign();
    end

    function setRho(uic, ~)
        if strcmp(uic.Style, 'slider')
            currentRho = round(uic.Value, 2);
            set(hRhoVal, 'String', currentRho)
        else
            currentRho = str2double(uic.String);
            set(hRho, 'Value', currentRho);
        end
        hClust.hCfg.log10RhoCut = currentRho;
        doReassign();
    end

    function setDelta(uic, ~)
        if strcmp(uic.Style, 'slider')
            currentDelta = round(uic.Value, 2);
            set(hDeltaVal, 'String', currentDelta)
        else
            currentDelta = str2double(uic.String);
            set(hDelta, 'Value', currentDelta);
        end
        hClust.hCfg.log10DeltaCut = currentDelta;
        doReassign();
    end

    function setSelectedClusters(uic, ~)
        selectedClusters = uic.Value;
        plotFeatures();
        plotWaves();
    end

    function doReassign()
        set(handles, 'Enable', 'off');
        hBox = msgbox('Reassigning clusters...please wait');
        hCFig.wait(1);
        hClust.reassign();
        updateSiteClusters();

        plotFeatures();
        plotRD();
        plotWaves();

        hCFig.wait(0);
        jrclust.utils.tryClose(hBox);
        set(handles, 'Enable', 'on');
    end

    function updateSiteClusters()
        % have to update the selected clusters
        siteSpikes = hClust.spikesBySite{currentSite};
        uniqueClusters = unique(hClust.spikeClusters(siteSpikes));
        uniqueClusters = uniqueClusters(uniqueClusters > 0);
        ucStr = arrayfun(@(iC) sprintf('Cluster %d', iC), uniqueClusters, 'UniformOutput', 0);

        set(hSiteClusters, 'Value', 1, 'String', '');

        selectedClusters = 1; % select all clusters

        set(hSiteClusters, 'String', ucStr, 'Min', 1, 'Max', 4);
        set(hSiteClusters, 'Value', selectedClusters);
    end

    function finalize(hObject, ~)
        hBox = msgbox('Reverting to your previous clustering...please wait');

        hClust.hCfg.log10RhoCut = initialSettings.rhoCut;
        hClust.hCfg.log10DeltaCut = initialSettings.deltaCut;
        hClust.hCfg.RDDetrendMode = initialSettings.RDDetrendMode;
        hClust.revert(initialSettings.position);

        hFigRD.close();
        hFigWav.close();
        jrclust.utils.tryClose(hBox);
        delete(hObject);
    end

	function plotFeatures()
        hCFig.cla();
        hCFig.axApply('default', @hold, 'on');

        %siteSpikes = [hClust.spikesBySite{site}; hClust.spikesBySite2{site}];
        siteSpikes = hClust.spikesBySite{currentSite};
        spikeClusters = hClust.spikeClusters(siteSpikes);

        cproj = strsplit(currentProjection, ' vs. ');
        % set XData
        if strcmp(cproj{2}, 'amp')
            XData = hClust.spikeAmps(siteSpikes);
            hCFig.axApply('default', @xlabel, 'Amplitude');
        elseif strcmp(cproj{2}, 'time')
            XData = hClust.spikeTimes(siteSpikes);
            hCFig.axApply('default', @xlabel, 'Time');
        elseif startsWith(cproj{2}, 'PC')
            nSitesRad = hClust.hCfg.nSitesEvt;
            whichPc = str2double(strrep(cproj{2}, 'PC', ''));
            iPc = (whichPc-1)*nSitesRad + whichPc;
            XData = squeeze(hClust.spikeFeatures(iPc, 1, hClust.spikesBySite{currentSite}));
            hCFig.axApply('default', @xlabel, cproj{2});
        end

        % set YData
        if strcmp(cproj{1}, 'amp')
            YData = hClust.spikeAmps(siteSpikes);
            hCFig.axApply('default', @ylabel, 'Amplitude');
        elseif startsWith(cproj{1}, 'PC')
            nSitesRad = hClust.hCfg.nSitesEvt;
            whichPc = str2double(strrep(cproj{1}, 'PC', ''));
            iPc = (whichPc-1)*nSitesRad + whichPc;
            YData = squeeze(hClust.spikeFeatures(iPc, 1, hClust.spikesBySite{currentSite}));
            hCFig.axApply('default', @ylabel, cproj{1});
        end

        uniqueClusters = cellfun(@(c) str2double(strrep(c, 'Cluster ', '')), ucStr, 'UniformOutput', 1);
        uniqueClusters = uniqueClusters(selectedClusters);
        for iCluster = 1:numel(uniqueClusters)
            cluster = uniqueClusters(iCluster);
            mask = jrclust.utils.subsample(find(spikeClusters == cluster), 1000); %#ok<*FNDSB>
            hCFig.addPlot(sprintf('s%dc%d', currentSite, cluster), XData(mask), YData(mask), '.');
        end
        hCFig.axApply('default', @hold, 'off');
        hCFig.toForeground();
        legend(ucStr(selectedClusters));
    end

    % plot rho-delta
    function plotRD()
        hFigRD.cla();

        if strcmp(hClust.hCfg.RDDetrendMode, 'global')
            [centers, rho, delta] = jrclust.cluster.densitypeaks.detrendRhoDelta(hClust, hClust.spikesBySite, 0, hClust.hCfg);
            delta = jrclust.utils.nanlog10(delta);
            fDetrend = 1;
        elseif strcmp(hClust.hCfg.RDDetrendMode, 'local')
            [centers, rho, delta] = jrclust.cluster.densitypeaks.detrendRhoDelta(hClust, hClust.spikesBySite, 1, hClust.hCfg);
            delta = jrclust.utils.nanlog10(delta);
            fDetrend = 1;
        else
            centers = find(hClust.spikeRho(:) > 10^(hClust.hCfg.log10RhoCut) & hClust.spikeDelta(:) > 10^(hClust.hCfg.log10DeltaCut));
            rho = jrclust.utils.nanlog10(hClust.spikeRho(:));
            delta = jrclust.utils.nanlog10(hClust.spikeDelta(:));
            fDetrend = 0;
        end

        mask = jrclust.utils.subsample(1:numel(rho), ceil(3*numel(rho)/5));

        hFigRD.addPlot('allSpikes', rho(mask), delta(mask), '.');
        hFigRD.axApply('default', @hold, 'on');

        hFigRD.axApply('default', @axis, 'tight');
        hFigRD.axApply('default', @axis, [min(jrclust.utils.nanlog10(hClust.spikeRho)) ...
                     max(jrclust.utils.nanlog10(hClust.spikeRho)) ...
                     min(jrclust.utils.nanlog10(hClust.spikeDelta)) ...
                     max(jrclust.utils.nanlog10(hClust.spikeDelta)) + 1]);
        hFigRD.axApply('default', @set, 'XScale', 'linear', 'YScale', 'linear');

        % show rho/delta cutoff lines
        hFigRD.addPlot('RDCuts', hClust.hCfg.log10RhoCut * [1 1], hFigRD.axApply('default', @get, 'YLim'), 'r--', ...
                       hFigRD.axApply('default', @get, 'XLim'), hClust.hCfg.log10DeltaCut*[1, 1], 'r--');
        hFigRD.axApply('default', @grid, 'on');

        % label cluster centers
        if ~isempty(hClust.clusterCenters)
            centers = hClust.clusterCenters; % do not overwrite
        end
        centersX = double(rho(centers));
        centersY = double(delta(centers));
        hFigRD.addPlot('centers', centersX, centersY, 'r.');

        % set labels
        hFigRD.axApply('default', @xlabel, 'log10 rho');
        if fDetrend
            hFigRD.axApply('default', @ylabel, 'log10 delta (detrended)');
        else
            hFigRD.axApply('default', @ylabel, 'log10 delta');
        end
        hFigRD.axApply('default', @title, sprintf('rho-cut: %f, delta-cut: %f', hClust.hCfg.log10RhoCut, hClust.hCfg.log10DeltaCut));

        drawnow;
    end

    function plotWaves()
        hFigWav.cla();
        hFigWav.axApply('default', @hold, 'on');

        siteSpikes = hClust.spikesBySite{currentSite};
        spikeClusters = hClust.spikeClusters(siteSpikes);

        uniqueClusters = cellfun(@(c) str2double(strrep(c, 'Cluster ', '')), ucStr, 'UniformOutput', 1);
        uniqueClusters = uniqueClusters(selectedClusters);

        YData = squeeze(hClust.spikesFilt(:, 1, siteSpikes));
        offset = max(max(YData) - min(YData)) + 100;
        nSamples = size(YData, 1);
        for iCluster = 1:numel(uniqueClusters)
            cluster = uniqueClusters(iCluster);
            mask = jrclust.utils.subsample(find(spikeClusters == cluster), 100);
            hFigWav.addPlot(sprintf('s%dc%d', currentSite, cluster), 1:nSamples, offset*(iCluster-1) + YData(:, mask));
        end
        hFigWav.axApply('default', @hold, 'off');
        hFigWav.toForeground();
        legend(ucStr(selectedClusters));
    end
end
