function updateSelect(obj, iClusters, force)
    %UPDATESELECT Select a (pair of) cluster(s) across all views
    iClusters = min(max(iClusters, 1), obj.hClust.nClusters);
    if numel(iClusters) > 2
        iClusters = iClusters(1:2);
    elseif isempty(iClusters)
        iClusters = 1;
    end
    if numel(iClusters) == 2 && diff(iClusters) == 0
        iClusters = iClusters(1);
    end

    % force update, e.g., after a delete
    if nargin < 3
        force = 0;
    end

    if ~force && jrclust.utils.isEqual(iClusters, obj.selected)
        return;
    end

    obj.selected = iClusters;

    iSite = obj.hClust.clusterSites(obj.selected(1));

    % update current site for amplitude view
    obj.currentSite = iSite;
    % limit the number of sites to display in the feature projection view
    nSites = min(obj.hCfg.nSitesFigProj, size(obj.hCfg.siteNeighbors, 1)); % by request

    % center sites around cluster center site
    if nSites < size(obj.hCfg.siteNeighbors, 1)
        obj.projSites = iSite:iSite + nSites - 1;
        if obj.projSites(end) > max(obj.hCfg.siteMap) % correct for overshooting
            obj.projSites = obj.projSites - max(obj.projSites) + max(obj.hCfg.siteMap);
        end
    else
        obj.projSites = sort(obj.hCfg.siteNeighbors(:, iSite), 'ascend');
    end
    
    obj.projSites = obj.hCfg.siteNeighbors(1:nSites,iSite);
    
    [~,idx] = sort(obj.channel_idx(obj.projSites));
    obj.projSites = obj.projSites(idx);

    % update plots
    obj.updateFigCorr();
    obj.updateFigHist();
    obj.updateFigISI();
    obj.updateFigMap();
    obj.updateFigPos();
    obj.updateFigProj(1);
    obj.updateFigTime(1);
    obj.updateFigPSTH(0);

    % update cursors
    obj.updateCursorFigWav();
    obj.updateCursorFigSim();

    % update menu
    obj.updateNoteMenu();
    obj.updateHistMenu();

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');

        % zoom if explicitly told to zoom
        doZoom = isfield(hFigWav.figData, 'zoom') && hFigWav.figData.zoom;

        % ... or if out of bounds
        if ~doZoom
            xRange = hFigWav.axApply('default', @xlim);
            yRange = hFigWav.axApply('default', @ylim);
            iCluster = obj.selected(1);
            if numel(obj.selected) > 1
                jCluster = obj.selected(2);
            else
                jCluster = iCluster;
            end

            iSite = obj.hClust.clusterSites(iCluster);
            jSite = obj.hClust.clusterSites(jCluster);

            % x values out of range
            doZoom = doZoom | obj.unitIndex(iCluster) < xRange(1) | obj.unitIndex(iCluster) > xRange(2);
            doZoom = doZoom | obj.unitIndex(jCluster) < xRange(1) | obj.unitIndex(jCluster) > xRange(2);

            % y values out of range
            doZoom = doZoom | iSite < yRange(1) | iSite > yRange(2);
            doZoom = doZoom | jSite < yRange(1) | jSite > yRange(2);
        end

        if doZoom
            obj.keyPressFigWav([], struct('Key', 'z'));
        end
    end
end