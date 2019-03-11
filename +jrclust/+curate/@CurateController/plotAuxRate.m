function plotAuxRate(obj, selectedOnly)
    %PLOTAUXRATE
    if selectedOnly
        corrData = jrclust.views.plotAuxCorr(obj.hClust, obj.selected(1));
    else
        corrData = jrclust.views.plotAuxCorr(obj.hClust, []);
    end

    if ~isempty(corrData)
        jrclust.utils.exportToWorkspace(corrData, obj.hCfg.verbose);
    end
end
