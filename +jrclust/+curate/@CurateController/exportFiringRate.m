function exportFiringRate(obj)
    %EXPORTFIRINGRATE
    firingRates = obj.hClust.firingRates();
    jrclust.utils.exportToWorkspace(struct('firingRates', firingRates), obj.hCfg.verbose);
end