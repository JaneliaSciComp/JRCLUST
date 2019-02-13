function saveBinaries(obj)
    %SAVEBINARIES Save raw/filtered traces, features to disk
    if isempty(obj.res)
        return;
    end
    sessionName = obj.hCfg.sessionName;

    % save spikesRaw
    if isfield(obj.res, 'spikesRaw') && ~isempty(obj.res.spikesRaw)
        filename = fullfile(obj.hCfg.outputDir, [sessionName, '_raw.jrc']);
        [flag, errmsg] = writeBin(filename, obj.res.spikesRaw, '*int16');
        if flag && obj.hCfg.verbose
            fprintf('Saved spikesRaw to %s\n', filename);
        elseif ~flag
            warning('Failed to save spikesRaw: %s', errmsg);
        end
    end

    % save spikesFilt
    if isfield(obj.res, 'spikesFilt') && ~isempty(obj.res.spikesFilt)
        filename = fullfile(obj.hCfg.outputDir, [sessionName, '_filt.jrc']);
        [flag, errmsg] = writeBin(filename, obj.res.spikesFilt, '*int16');
        if flag && obj.hCfg.verbose
            fprintf('Saved spikesFilt to %s\n', filename);
        elseif ~flag
            warning('Failed to save spikesFilt: %s', errmsg);
        end
    end

    % save spikeFeatures
    if isfield(obj.res, 'spikeFeatures') && ~isempty(obj.res.spikeFeatures)
        filename = fullfile(obj.hCfg.outputDir, [sessionName, '_features.jrc']);
        [flag, errmsg] = writeBin(filename, obj.res.spikeFeatures, '*single');
        if flag && obj.hCfg.verbose
            fprintf('Saved spikeFeatures to %s\n', filename);
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