function loadFiles(obj)
    %LOADFILES Load results struct
    if obj.isError
        error(obj.errMsg);
    end

    if ~exist(obj.hCfg.resFile, 'file')
        warning('%s does not exist', obj.hCfg.resFile);
        return;
    end

    try
        if obj.hCfg.verbose
            fprintf('Loading %s...', obj.hCfg.resFile);
            t = tic;
        end
        res_ = load(obj.hCfg.resFile);
        if obj.hCfg.verbose
            fprintf('done (took %0.2f s)\n', toc(t));
        end
    catch ME
        warning('failed to load %s: %s', ME.message);
        return;
    end

    if isfield(res_, 'spikeTimes')
        % load spikesRaw
        if isfield(res_, 'rawShape')
            if obj.hCfg.verbose
                t = tic;
                fprintf('Loading %s...', obj.hCfg.rawFile);
            end
            spikesRaw = readBin(obj.hCfg.rawFile, res_.rawShape, '*int16');
            if obj.hCfg.verbose
                fprintf('done (took %0.2f s)\n', toc(t));
            end
        else
            spikesRaw = [];
        end

        % load spikesFilt
        if isfield(res_, 'filtShape')
            if obj.hCfg.verbose
                t = tic;
                fprintf('Loading %s...', obj.hCfg.filtFile);
            end
            spikesFilt = readBin(obj.hCfg.filtFile, res_.filtShape, '*int16');
            if obj.hCfg.verbose
                fprintf('done (took %0.2f s)\n', toc(t));
            end
        else
            spikesFilt = [];
        end

        if isfield(res_, 'featuresShape')
            if obj.hCfg.verbose
                t = tic;
                fprintf('Loading %s...', obj.hCfg.featuresFile);
            end
            spikeFeatures = readBin(obj.hCfg.featuresFile, res_.featuresShape, '*single');
            if obj.hCfg.verbose
                fprintf('done (took %0.2f s)\n', toc(t));
            end
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