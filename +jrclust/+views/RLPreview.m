function hCFig = RLPreview(hClust)
    hCFig = jrclust.views.ControlFigure('FigRLPreview', [0.1 0.1 0.5 0.5], 'RL Preview', true, false);
    hCFig.axes();

    % populate the site # dropdown
    hCFig.addUicontrol('siteNoLabel', 'Style', 'text', ...
                       'String', 'Site #', ...
                       'Position', [10, 0, 30, 25]);
    hCFig.addUicontrol('siteNo', 'Style', 'popup', ...
                       'String', arrayfun(@(n) num2str(n), 1:numel(hClust.spikesBySite), 'UniformOutput', false), ...
                       'Position', [45, 0, 50, 30], ...
                       'Callback', @setSite);

    % populate the feature projection dropdown
    if strcmp(hClust.hCfg.clusterFeature, 'pca')
        nPCs = hClust.hCfg.nPcPerChan;
        projections = cell(nPCs, 1);
        for i = 1:nPCs
            projections{i} = sprintf('PC%d vs. amp', i);
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
                       'Position', [125, 0, 50, 25]);
    hCFig.addUicontrol('projection', 'Style', 'popup', ...
                       'String', projections, ...
                       'Position', [175, 0, 100, 30], ...
                       'Callback', @setProjection);
    hCFig.title('Projections of spikes onto features');
    hFigRD = jrclust.views.Figure('FigRD', [0.6 0.1 0.4 0.5], ['Cluster rho-delta: ', hClust.hCfg.sessionName], false, false);

    % populate the detrend option dropdown
    detrendOptions = {'none', 'local', 'global', 'hidehiko'};
    others = detrendOptions(~strcmp(detrendOptions, hClust.hCfg.rlDetrendMode));
    detrendOptions{1} = hClust.hCfg.rlDetrendMode;
    detrendOptions(2:end) = others;

    hCFig.addUicontrol('detrendLabel', 'Style', 'text', ...
                       'String', 'Projection', ...
                       'Position', [285, 0, 50, 25]);
    hCFig.addUicontrol('projection', 'Style', 'popup', ...
                       'String', detrendOptions, ...
                       'Position', [345, 0, 100, 30], ...
                       'Callback', @setDetrend);

    % populate the rho slider
    rhoMin = min(hClust.spikeRho);
    rhoMax = max(hClust.spikeRho);

%     hCFig.addUicontrol('detrendLabel', 'Style', 'text', ...
%                        'String', 'Projection', ...
%                        'Position', [275, 0, 50, 25]);
%     hCFig.addUicontrol('projection', 'Style', 'popup', ...
%                        'String', detrendOptions, ...
%                        'Position', [325, 0, 100, 30], ...
%                        'Callback', @setDetrend);

    % initial state
    initialSettings = struct('rlDetrendMode', hClust.hCfg.rlDetrendMode, ...
                             'position', size(hClust.history, 1) - 1, ...
                             'rhoCut', hClust.hCfg.log10RhoCut, ...
                             'deltaCut', hClust.hCfg.log10DeltaCut);
    currentSite = 1;
    currentProjection = projections{1};
    currentDetrend = detrendOptions{1};

    plotFeatures();
    plotRD();

    function setSite(uic, ~)
        if currentSite == uic.Value
            return;
        end
        currentSite = uic.Value;
        plotFeatures();
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
        hClust.hCfg.rlDetrendMode = currentDetrend;
        hCFig.wait(true);
        hClust.reassign();
        plotFeatures();
        plotRD();
        hCFig.wait(false);
    end

    function finalize(hObject, ~)
        hClust.revert(initialSettings.position);
        hClust.hCfg.log10RhoCut = initialSettings.rhoCut;
        hClust.hCfg.log10DeltaCut = initialSettings.deltaCut;
        hClust.hCfg.rlDetrendMode = initialSettings.rlDetrendMode;

        hFigRD.close();
        delete(hObject);
    end

	function plotFeatures()
        hCFig.cla();
        hCFig.hold('on');

        %siteSpikes = [hClust.spikesBySite{site}; hClust.spikesBySite2{site}];
        siteSpikes = hClust.spikesBySite{currentSite};
        spikeClusters = hClust.spikeClusters(siteSpikes);
        %spikeFeatures = [squeeze(hClust.spikeFeatures(1, 1, hClust.spikesBySite{site})); squeeze(hClust.spikeFeatures(1, 2, hClust.spikesBySite2{site}))];
        cproj = strsplit(currentProjection, ' vs. ');
        if strcmp(cproj{2}, 'amp')
            XData = hClust.spikeAmps(siteSpikes);
            hCFig.xlabel('Amplitude');
        elseif startsWith(cproj{2}, 'PC')
            nSitesRad = hClust.hCfg.nSitesEvt;
            whichPc = str2double(strrep(cproj{2}, 'PC', ''));
            iPc = (whichPc-1)*nSitesRad + whichPc;
            XData = squeeze(hClust.spikeFeatures(iPc, 1, hClust.spikesBySite{currentSite}));
            hCFig.xlabel(cproj{2});
        end

        if startsWith(cproj{1}, 'PC')
            nSitesRad = hClust.hCfg.nSitesEvt;
            whichPc = str2double(strrep(cproj{1}, 'PC', ''));
            iPc = (whichPc-1)*nSitesRad + whichPc;
            YData = squeeze(hClust.spikeFeatures(iPc, 1, hClust.spikesBySite{currentSite}));
            hCFig.ylabel(cproj{1});
        end

        uniqueClusters = unique(spikeClusters);
        uniqueClusters = uniqueClusters(uniqueClusters > 0);
        for iCluster = 1:numel(uniqueClusters)
            cluster = uniqueClusters(iCluster);
            mask = (spikeClusters == cluster);
            hCFig.addPlot(sprintf('s%dc%d', currentSite, cluster), XData(mask), YData(mask), '.');
        end
        hCFig.hold('off');
        hCFig.toForeground();
        set(gcf, 'CloseRequestFcn', @finalize);
        legendNames = arrayfun(@(iC) sprintf('Cluster %d', iC), uniqueClusters, 'UniformOutput', false);
        legend(legendNames);
    end

    % plot rho-delta
    function plotRD()
        hFigRD.cla();
        if strcmp(hClust.hCfg.rlDetrendMode, 'global')
            [centers, rho, delta] = jrclust.cluster.densitypeaks.detrendRhoDelta(hClust, hClust.spikesBySite, false, hClust.hCfg);
            delta = jrclust.utils.nanlog10(delta);
            fDetrend = true;
        elseif strcmp(hClust.hCfg.rlDetrendMode, 'local')
            [centers, rho, delta] = jrclust.cluster.densitypeaks.detrendRhoDelta(hClust, hClust.spikesBySite, true, hClust.hCfg);
            delta = jrclust.utils.nanlog10(delta);
            fDetrend = true;
        else
            centers = find(hClust.spikeRho(:) > 10^(hClust.hCfg.log10RhoCut) & hClust.spikeDelta(:) > 10^(hClust.hCfg.log10DeltaCut));
            rho = jrclust.utils.nanlog10(hClust.spikeRho(:));
            delta = jrclust.utils.nanlog10(hClust.spikeDelta(:));
            fDetrend = false;
        end

        hFigRD.addPlot('allSpikes', rho, delta, '.');
        hFigRD.hold('on');

        hFigRD.axis('tight');
        hFigRD.axis([-4 -.5 -1 2]);
        hFigRD.axSet('XScale', 'linear', 'YScale', 'linear');

        % show rho/delta cutoff lines
        hFigRD.addPlot('RDCuts', hClust.hCfg.log10RhoCut * [1 1], hFigRD.axGet('YLim'), 'r--', ...
                       hFigRD.axGet('XLim'), hClust.hCfg.log10DeltaCut*[1, 1], 'r--');
        hFigRD.grid('on');

        % label cluster centers
        if ~isempty(hClust.clusterCenters)
            centers = hClust.clusterCenters; % do not overwrite
        end
        centersX = double(rho(centers));
        centersY = double(delta(centers));
        hFigRD.addPlot('centers', centersX, centersY, 'r.');

        % set labels
        hFigRD.xlabel('log10 rho');
        if fDetrend
            hFigRD.ylabel('log10 delta (detrended)');
        else
            hFigRD.ylabel('log10 delta');
        end
        hFigRD.title(sprintf('rho-cut: %f, delta-cut: %f', hClust.hCfg.log10RhoCut, hClust.hCfg.log10DeltaCut));

        drawnow;
    end
end