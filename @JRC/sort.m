function sRes = sort(obj)
    %SORT Cluster detected spikes
    if isempty(obj.res)
        obj.loadFiles();
    end
    obj.isSort = 1;

    if ~isfield(obj.res, 'spikeFeatures')
        dlgAns = questdlg('Could not find all required data. Detect?', 'Detection required', 'No');
        if strcmp(dlgAns, 'Yes')
            obj.detect();
        else
            obj.isCompleted = 1;
            return;
        end
    end

    if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
        obj.hCfg.updateLog('detectedOn', sprintf('Using spikes detected on %s', datestr(obj.res.detectedOn)), 0, 0);
    elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
        obj.hCfg.updateLog('detectedOn', 'Using previously-detected spikes', 0, 0);
    end

    % clear GPU memory and set random seeds
    obj.clearMemory();

    % start the parallel pool
    if obj.hCfg.useParfor
        obj.startParPool();
    end

    obj.hCfg.updateLog('sortStep', 'Sorting detected spikes', 1, 0);
    obj.hSort = jrclust.sort.SortController(obj.hCfg);
    sRes = obj.hSort.sort(obj.res);

    if obj.hSort.isError
        obj.error(obj.hSort.errMsg);
    end
    obj.hCfg.updateLog('sortStep', 'Finished sorting', 0, 1);

    obj.res = jrclust.utils.mergeStructs(obj.res, sRes);
    obj.saveRes(obj.isDetect); % force overwrite if we're detecting

    obj.summarize();
end