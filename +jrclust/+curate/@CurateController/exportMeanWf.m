function exportMeanWf(obj, exportAll)
    %EXPORTMEANWF Export mean waveforms to workspace
    if exportAll
        jrclust.utils.exportToWorkspace(struct('meanWfGlobal', obj.hClust.meanWfGlobal), obj.hCfg.verbose);
    elseif ~isempty(obj.selected)
        mwf = obj.hClust.meanWfGlobal;

        primarySites = obj.hCfg.siteNeighbors(:, obj.hClust.clusterSites(obj.selected(1)));
        selectedWf = mwf(:, primarySites, obj.selected(1));
        if numel(obj.selected) == 2
            primarySites2 = obj.hCfg.siteNeighbors(:, obj.hClust.clusterSites(obj.selected(2)));
            selectedWf2 = mwf(:, primarySites2, obj.selected(2));
            jrclust.utils.exportToWorkspace(struct('selectedWf', selectedWf, 'selectedWf2', selectedWf2), obj.hCfg.verbose);
        else
            jrclust.utils.exportToWorkspace(struct('selectedWf', selectedWf), obj.hCfg.verbose);
        end
    end
end