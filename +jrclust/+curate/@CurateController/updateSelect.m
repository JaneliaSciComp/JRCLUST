function updateSelect(obj, iClusters)
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
    obj.updateMenu();
end