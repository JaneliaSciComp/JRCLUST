function clearMemory(obj)
    %CLEARMEMORY Clear GPU memory and set random seeds
    if obj.hCfg.useGPU
        % while we're here, clear GPU memory
        if obj.isDetect || obj.isSort
            if obj.hCfg.verbose
                fprintf('Clearing GPU memory...');
            end
            gpuDevice(); % selects GPU device
            gpuDevice([]); % clears GPU memory
            if obj.hCfg.verbose
                fprintf('done\n');
            end
        end

        parallel.gpu.rng(obj.hCfg.randomSeed);
    end

    rng(obj.hCfg.randomSeed);
end