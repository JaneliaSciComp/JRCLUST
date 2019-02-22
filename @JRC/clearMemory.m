function clearMemory(obj)
    %CLEARMEMORY Clear GPU memory and set random seeds
    if obj.hCfg.useGPU
        if obj.isDetect || obj.isSort
            obj.hCfg.updateLog('gpuMemory', 'Clearing GPU memory', 1, 0);
            gpuDevice(); % selects GPU device
            gpuDevice([]); % clears GPU memory
            obj.hCfg.updateLog('gpuMemory', 'GPU memory cleared', 0, 1);
        end

        parallel.gpu.rng(obj.hCfg.randomSeed);
    end

    rng(obj.hCfg.randomSeed);
end