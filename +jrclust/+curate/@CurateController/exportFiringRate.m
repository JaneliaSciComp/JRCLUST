function exportFiringRate(obj)
    %EXPORTFIRINGRATE
    firingRates = obj.hClust.getFiringRates);
    jrclust.utils.exportToWorkspace(struct('firingRates', firingRates), obj.hCfg.verbose);
end