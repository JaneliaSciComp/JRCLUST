function res = doLoadFiles(hCfg)
    %DOLOADFILES Load results struct
    res = struct();

    filename = fullfile(hCfg.outputDir, [hCfg.sessionName '_res.mat']);
    if ~isfile(filename)
        warning('%s does not exist', filename);
        return;
    end

    try
        res = load(filename);
    catch ME
        warning(ME.identifier, 'failed to load %s: %s', ME.message);
        return;
    end

    if isfield(res, 'spikeTimes')
        % load spikesRaw
        if isfield(res, 'rawShape')
            rawFile = fullfile(hCfg.outputDir, [hCfg.sessionName '_raw', '.jrc']);
            spikesRaw = loadBin(rawFile, res.rawShape, 'int16');
        else
            spikesRaw = [];
        end

        % load spikesFilt
        if isfield(res, 'filtShape')
            filtFile = fullfile(hCfg.outputDir, [hCfg.sessionName '_filt', '.jrc']);
            spikesFilt = loadBin(filtFile, res.filtShape, 'int16');
        else
            spikesFilt = [];
        end

        if isfield(res, 'featuresShape')
            featuresFile = fullfile(hCfg.outputDir, [hCfg.sessionName '_features', '.jrc']);
            spikeFeatures = loadBin(featuresFile, res.featuresShape, 'single');
        else
            spikeFeatures = [];
        end

        if isfield(res, 'hClust')
            res.hClust.spikesRaw = spikesRaw;
            res.hClust.spikesFilt = spikesFilt;
            res.hClust.spikeFeatures = spikeFeatures;
        else
            res.spikesRaw = spikesRaw;
            res.spikesFilt = spikesFilt;
            res.spikeFeatures = spikeFeatures;
        end
    end
end

%% LOCAL FUNCTIONS
function binData = loadBin(filename, binShape, dataType)
    %LOADBIN Save traces/features to binary file
    if isfile(filename)
        fid = fopen(filename, 'r');
        binData = fread(fid, Inf, ['*' dataType]);
        fclose(fid);
        binData = reshape(binData, binShape);
    else
        binData = [];
    end
end