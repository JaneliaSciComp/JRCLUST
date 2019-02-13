function loadFiles(obj)
    %LOADFILES Load results struct
    if obj.isError
        error(obj.errMsg);
    end

    filename = fullfile(obj.hCfg.outputDir, [obj.hCfg.sessionName '_res.mat']);
    if ~exist(filename, 'file')
        warning('%s does not exist', filename);
        return;
    end

    try
        res_ = load(filename);
    catch ME
        warning('failed to load %s: %s', ME.message);
        return;
    end

    if isfield(res_, 'spikeTimes')
        % load spikesRaw
        if isfield(res_, 'rawShape')
            rawFile = fullfile(obj.hCfg.outputDir, [obj.hCfg.sessionName '_raw', '.jrc']);
            spikesRaw = readBin(rawFile, res_.rawShape, '*int16');
        else
            spikesRaw = [];
        end

        % load spikesFilt
        if isfield(res_, 'filtShape')
            filtFile = fullfile(obj.hCfg.outputDir, [obj.hCfg.sessionName '_filt', '.jrc']);
            spikesFilt = readBin(filtFile, res_.filtShape, '*int16');
        else
            spikesFilt = [];
        end

        if isfield(res_, 'featuresShape')
            featuresFile = fullfile(obj.hCfg.outputDir, [obj.hCfg.sessionName '_features', '.jrc']);
            spikeFeatures = readBin(featuresFile, res_.featuresShape, '*single');
        else
            spikeFeatures = [];
        end

        % set spikesRaw/spikesFilt/spikeFeatures (warn if empty!)
        if isempty(spikesRaw)
            warning('spikesRaw is empty');
        end
        res_.spikesRaw = spikesRaw;

        if isempty(spikesFilt)
            warning('spikesFilt is empty');
        end
        res_.spikesFilt = spikesFilt;

        if isempty(spikeFeatures)
            warning('spikeFeatures is empty');
        end
        res_.spikeFeatures = spikeFeatures;

        if isfield(res_, 'hClust')
            res_.hClust.spikesRaw = spikesRaw;
            res_.hClust.spikesFilt = spikesFilt;
            res_.hClust.spikeFeatures = spikeFeatures;
        end
    end

    obj.res = res_;
    if isfield(obj.res, 'hClust')
        obj.res.hClust.hCfg = obj.hCfg;
    end
end

%% LOCAL FUNCTIONS
function binData = readBin(filename, binShape, dataType)
    %LOADBIN Save traces/features to binary file
    if exist(filename, 'file')
        fid = fopen(filename, 'r');
        binData = fread(fid, Inf, dataType);
        fclose(fid);
        binData = reshape(binData, binShape);
    else
        binData = [];
    end
end