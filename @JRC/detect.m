function dRes = detect(obj)
    %DETECT Detect spikes in recording
    obj.res = []; % reset obj.res
    obj.isDetect = 1;

    % clear GPU memory and set random seeds
    obj.clearMemory();

    % start the parallel pool
    if obj.hCfg.useParfor
        obj.startParPool();
    end

    obj.hDetect = jrclust.detect.DetectController(obj.hCfg);
    dRes = obj.hDetect.detect();

    if obj.hDetect.isError
        obj.error(obj.hDetect.errMsg);
    elseif obj.hCfg.verbose
        fprintf('Detection completed in %0.2f seconds\n', dRes.detectTime);
    end

    obj.res = dRes;

    % save files
    obj.saveBinaries();
    obj.saveRes(0);
end