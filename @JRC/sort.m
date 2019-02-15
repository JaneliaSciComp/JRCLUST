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

    if obj.hCfg.verbose
        % inform user we're using previously detected spikes
        if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
            fprintf('Using spikes detected on %s\n', datestr(obj.res.detectedOn));
        elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
            fprintf('Using previously-detected spikes\n');
        end
    end

    % clear GPU memory and set random seeds
    obj.clearMemory();

    % start the parallel pool
    if obj.hCfg.useParfor
        obj.startParPool();
    end

    obj.hSort = jrclust.sort.SortController(obj.hCfg);
    sRes = obj.hSort.sort(obj.res);

    if obj.hSort.isError
        obj.error(obj.hSort.errMsg);
    elseif obj.hCfg.verbose
        fprintf('Sorting completed in %0.2f seconds\n', sRes.sortTime);
    end

    obj.res = jrclust.utils.mergeStructs(obj.res, sRes);
    obj.saveRes(obj.isDetect); % force overwrite if we're detecting
end