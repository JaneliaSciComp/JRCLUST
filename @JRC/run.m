function run(obj)
    %RUN Run commands
    if obj.isError
        error(obj.errMsg);  
    elseif obj.isCompleted
        warning('command ''%s'' completed successfully; to rerun, use rerun()', obj.cmd);
        return;
    end

    % DETECT SPIKES
    gpuDetect = obj.hCfg.useGPU; % save this in case useGPU is disabled during detection step
    if obj.isDetect
        obj.detect();
    end

    % CLUSTER SPIKES
    gpuSort = obj.hCfg.useGPU | gpuDetect;
    if obj.isSort
        obj.hCfg.useGPU = gpuSort;
        obj.sort();
    end

    % CURATE SPIKES
    gpuCurate = obj.hCfg.useGPU | gpuSort;
    if obj.isCurate
        obj.hCfg.useGPU = gpuCurate;
        obj.curate();
    end

    obj.hCfg.closeLog();
    obj.isCompleted = 1;

    % save full params for reproducibility
    if ~obj.hCfg.getOr('testRun', 0)
        [~, ~, ext] = fileparts(obj.hCfg.configFile);
        newExt = sprintf('-full%s', ext);
        fullFilename = jrclust.utils.subsExt(obj.hCfg.configFile, newExt);
        obj.hCfg.save(fullFilename, 1);
    end
end