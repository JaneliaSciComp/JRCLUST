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

    obj.hCfg.updateLog('detectStep', 'Detecting spikes in recordings', 1, 0);
    obj.hDetect = jrclust.detect.DetectController(obj.hCfg);
    dRes = obj.hDetect.detect();

    if obj.hDetect.isError
        obj.error(obj.hDetect.errMsg);
    end
    obj.hCfg.updateLog('detectStep', 'Finished detecting', 0, 1);

    obj.res = dRes;

    % save files
    obj.saveRes(0);

    if ~obj.isSort
        obj.summarize();
    end
end
