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

    if obj.hCfg.verbose
        % inform user we're using previously detected spikes
        if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
            fprintf('Using spikes detected on %s\n', datestr(obj.res.detectedOn));
        elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
            fprintf('Using previously-detected spikes\n');
        end

        % inform user we're using a previously-computed clustering
        if ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'sortedOn')
            fprintf('Using clustering computed on %s\n', datestr(obj.res.sortedOn));
        elseif ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'hClust')
            fprintf('Using previously-clustered spikes\n');
        end

        % inform user of the last time this set was curated
        if ~isempty(obj.res) && isfield(obj.res, 'curatedOn')
            fprintf('Last manually edited on %s\n', datestr(obj.res.curatedOn));
        end
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