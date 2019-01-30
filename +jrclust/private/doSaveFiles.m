function doSaveFiles(res, hCfg)
    %DOSAVEFILES Save results structs and binary files to disk
    if isempty(res)
        return;
    end

    filename = fullfile(hCfg.outputDir, [hCfg.sessionName '_res.mat']);
    if exist(filename, 'file') % don't overwrite without explicit permission
        question = sprintf('%s already exists. Overwrite?', filename);
        dlgAns = questdlg(question, 'Confirm overwrite', 'No');
        if isempty(dlgAns) || ismember(dlgAns, {'No', 'Cancel'})
            return;
        end
    end

    % save spikesRaw
    res = saveBin(res, 'spikesRaw', 'raw', hCfg);
    if isfield(res, 'spikesRaw') % save to file failed
        warning('Failed to save spikesRaw; will be saved in %s', filename);
    end

    % save spikesFilt
    res = saveBin(res, 'spikesFilt', 'filt', hCfg);
    if isfield(res, 'spikesFilt') % save to file failed
        warning('Failed to save spikesFilt; will be saved in %s', filename);
    end

    % save spikeFeatures
    res = saveBin(res, 'spikeFeatures', 'features', hCfg);
    if isfield(res, 'spikeFeatures') % save to file failed
        warning('Failed to save spikeFeatures; will be saved in %s', filename);
    end

    % save everything else
    if isfield(res, 'hClust') && ~isempty(res.hClust.spikesRaw)
        doRestore = 1;
        spikesRaw = res.hClust.spikesRaw;
        res.hClust.spikesRaw = [];
        spikesFilt = res.hClust.spikesFilt;
        res.hClust.spikesFilt = [];
        spikeFeatures = res.hClust.spikeFeatures;
        res.hClust.spikeFeatures = [];
    else
        doRestore = 0;
    end

    jrclust.utils.saveStruct(res, filename);
    if hCfg.verbose
        fprintf('Saved results to %s\n', filename);
    end

    if doRestore
        res.hClust.spikesRaw = spikesRaw;
        res.hClust.spikesFilt = spikesFilt;
        res.hClust.spikeFeatures = spikeFeatures;
    end
end

%% LOCAL FUNCTIONS
function res = saveBin(res, binField, tag, hCfg)
    %SAVEBIN Save traces/features to binary file
    if isfield(res, binField)
        binFile = fullfile(hCfg.outputDir, [hCfg.sessionName '_', tag, '.jrc']);
        fid = fopen(binFile, 'w');
        if fid ~= -1
            prec = jrclust.utils.ifEq(strcmp(tag, 'features'), '*single', '*int16');
            ct = fwrite(fid, res.(binField), prec);
            fclose(fid);

            if ct ~= numel(res.(binField)) % write failed
                return;
            end

            if hCfg.verbose
                fprintf('Saved %s to %s\n', binField, binFile);
            end

            res.([tag 'Shape']) = size(res.(binField));
            res = rmfield(res, binField);
        end
    end
end