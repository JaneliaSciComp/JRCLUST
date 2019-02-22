function curate(obj)
    %CURATE Spin up the manual GUI for curation
    if isempty(obj.res)
        obj.loadFiles();
    end
    obj.isCurate = 1;

    if ~isfield(obj.res, 'hClust')
        dlgAns = questdlg('Could not find all required data. Sort?', 'Sorting required', 'No');
        if strcmp(dlgAns, 'Yes')
            obj.sort();
        else
            obj.isCompleted = 1;
            return;
        end
    end

    % inform user we're using previously detected spikes
    if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
        obj.hCfg.updateLog('detectedOn', sprintf('Using spikes detected on %s', datestr(obj.res.detectedOn)), 0, 0);
    elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
        obj.hCfg.updateLog('detectedOn', 'Using previously-detected spikes');
    end

    % inform user we're using a previously-computed clustering
    if ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'sortedOn')
        obj.hCfg.updateLog('sortedOn', sprintf('Using clustering computed on %s', datestr(obj.res.sortedOn)), 0, 0);
    elseif ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'hClust')
        obj.hCfg.updateLog('sortedOn', 'Using previously-clustered spikes', 0, 0);
    end

    % inform user of the last time this set was curated
    if ~isempty(obj.res) && isfield(obj.res, 'curatedOn')
        obj.hCfg.updateLog('curatedOn', sprintf('Last manually edited on %s\n', datestr(obj.res.curatedOn)), 0, 0);
    end

    % clear GPU memory and set random seeds
    obj.clearMemory();

    % start the parallel pool
    if obj.hCfg.useParfor
        obj.startParPool();
    end

    obj.hCurate = jrclust.curate.CurateController(obj.res);
    obj.hCurate.beginSession();
end