function beginSession(obj)
    %BEGINSESSION Start curating clusters
    if ~isempty(obj.hFigs) % session already running
        return;
    end

    obj.isEnding = 0;
    obj.cRes = struct('hClust', obj.hClust);
    obj.maxAmp = obj.hCfg.maxAmp;

    % select first cluster
    obj.selected = 1;
    obj.currentSite = obj.hClust.clusterSites(1);

    obj.plotAllFigures();
end