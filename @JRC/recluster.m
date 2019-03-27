function recluster(obj)
    %RECLUSTER Recluster spikes
    if isempty(obj.res)
        obj.loadFiles();
    end

    if isfield(obj.res, 'hClust')
        obj.isSort = 1;
        obj.res.hClust.hCfg = obj.hCfg; % update hClust's config
        obj.res.hClust.reassign();
        obj.res.hClust.autoMerge();
        obj.res.sortedOn = now();
        obj.saveRes(1);

        obj.summarize();
    else
        obj.isError = 1;
        obj.errMsg = 'hClust not found';
    end
end