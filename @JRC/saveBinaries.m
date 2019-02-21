function saveBinaries(obj)
    %SAVEBINARIES Save raw/filtered traces, features to disk
    if isempty(obj.res)
        return;
    end

    % save spikesRaw
    if isfield(obj.res, 'spikesRaw') && ~isempty(obj.res.spikesRaw)
        obj.hCfg.updateLog('saveBin', sprintf('Saving spikesRaw to %s', obj.hCfg.rawFile), 1, 0);
        [flag, errmsg] = writeBin(obj.hCfg.rawFile, obj.res.spikesRaw, '*int16');
        if flag
            obj.hCfg.updateLog('saveBin', 'Finished saving spikesRaw', 0, 1);
        elseif ~flag
            warning('Failed to save spikesRaw: %s', errmsg);
        end
    end

    % save spikesFilt
    if isfield(obj.res, 'spikesFilt') && ~isempty(obj.res.spikesFilt)
        obj.hCfg.updateLog('saveBin', sprintf('Saving spikesFilt to %s', obj.hCfg.filtFile), 1, 0);
        [flag, errmsg] = writeBin(obj.hCfg.filtFile, obj.res.spikesFilt, '*int16');
        if flag
            obj.hCfg.updateLog('saveBin', 'Finished saving spikesFilt', 0, 1);
        elseif ~flag
            warning('Failed to save spikesFilt: %s', errmsg);
        end
    end

    % save spikeFeatures
    if isfield(obj.res, 'spikeFeatures') && ~isempty(obj.res.spikeFeatures)
        obj.hCfg.updateLog('saveBin', sprintf('Saving spikeFeatures to %s', obj.hCfg.featuresFile), 1, 0);
        [flag, errmsg] = writeBin(obj.hCfg.featuresFile, obj.res.spikeFeatures, '*single');
        if flag
            obj.hCfg.updateLog('saveBin', 'Finished saving spikeFeatures', 0, 1);
        elseif ~flag
            warning('Failed to save spikeFeatures: %s', errmsg);
        end
    end
end

%% LOCAL FUNCTIONS
function [success, msg] = writeBin(filename, data, dataType)
    [fid, msg] = fopen(filename, 'w');
    if fid == -1
        success = 0;
    else
        ct = fwrite(fid, data, dataType);
        fclose(fid);
        success = (ct == numel(data));
        msg = sprintf('%d/%d elements written', ct, numel(data));
    end
end