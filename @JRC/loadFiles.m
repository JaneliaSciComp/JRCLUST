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
        obj.hCfg.updateLog('loadRes', sprintf('Loading %s', obj.hCfg.resFile), 1, 0);
        res_ = load(obj.hCfg.resFile);
        obj.hCfg.updateLog('loadRes', sprintf('Finished loading %s', obj.hCfg.resFile), 0, 1);
    catch ME
        warning('Failed to load %s: %s', ME.message);
        return;
    end

    if isfield(res_, 'spikeTimes')
        % load spikesRaw
        if isfield(res_, 'rawShape')
            obj.hCfg.updateLog('loadRaw', sprintf('Loading %s', obj.hCfg.rawFile), 1, 0);
            spikesRaw = readBin(obj.hCfg.rawFile, res_.rawShape, '*int16');
            obj.hCfg.updateLog('loadRaw', sprintf('Finished loading %s', obj.hCfg.rawFile), 0, 1);
        else
            spikesRaw = [];
        end

        % load spikesFilt
        if isfield(res_, 'filtShape')
            obj.hCfg.updateLog('loadFilt', sprintf('Loading %s', obj.hCfg.filtFile), 1, 0);
            spikesFilt = readBin(obj.hCfg.filtFile, res_.filtShape, '*int16');
            obj.hCfg.updateLog('loadFilt', sprintf('Finished loading %s', obj.hCfg.filtFile), 0, 1);
        else
            spikesFilt = [];
        end

        % load spikeFeatures
        if isfield(res_, 'featuresShape')
            obj.hCfg.updateLog('loadFeatures', sprintf('Loading %s', obj.hCfg.featuresFile), 1, 0);
            spikeFeatures = readBin(obj.hCfg.featuresFile, res_.featuresShape, '*single');
            obj.hCfg.updateLog('loadFeatures', sprintf('Finished loading %s', obj.hCfg.featuresFile), 0, 1);
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

        % restore values to hClust
        if isfield(res_, 'hClust')
            if isa(res_.hClust, 'jrclust.models.clustering.DensityPeakClustering')
                res_.hClust = convertToNew(res_, obj.hCfg);
            end

            hClustFields = fieldnames(res_.hClust);
            for i = 1:numel(hClustFields)
                fn = hClustFields{i};
                if isempty(res_.hClust.(fn)) && isfield(res_, fn)
                    res_.hClust.(fn) = res_.(fn);
                elseif ismember(fn, {'clusterCenters', 'clusterCentroids'}) && isfield(res_, fn)
                    res_.hClust.sRes.(fn) = res_.(fn);
                end
            end

            % restore initialClustering
            res_.hClust.initialClustering = res_.spikeClusters;

            % supply hClust with our own hCfg
            res_.hClust.hCfg = obj.hCfg;
        elseif isfield(res_, 'spikeTemplates') % create a new TemplateClustering
            hClust = jrclust.sort.TemplateClustering(struct(), struct(), obj.hCfg);
            fieldNames = fieldnames(res_);
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                if isprop(hClust, fn)
                    hClust.(fn) = res_.(fn);
                end
            end

            res_.hClust = hClust;
        elseif isfield(res_, 'spikeClusters')
            hClust = jrclust.sort.DensityPeakClustering(struct(), struct(), obj.hCfg);
            fieldNames = fieldnames(res_);
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                if isprop(hClust, fn)
                    hClust.(fn) = res_.(fn);
                end
            end

            res_.hClust = hClust;
        end

        if isfield(res_, 'hRecs') % don't try to load recordings
            res_ = rmfield(res_, 'hRecs');
        end
    else
        warning('spikeTimes not found in %s', obj.hCfg.resFile);
        res_ = struct();
    end

    obj.res = res_;
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

function hClustNew = convertToNew(res, hCfg)
    hClustOld = res.hClust;

    dRes = hClustOld.dRes;
    dRes.spikesRaw = res.spikesRaw;
    dRes.spikesFilt = res.spikesFilt;
    dRes.spikeFeatures = res.spikeFeatures;

    sRes = hClustOld.sRes;
    if isfield(sRes, 'simScore')
        sRes.waveformSim = sRes.simScore;
        sRes = rmfield(sRes, 'simScore');
    end

    hClustNew = jrclust.sort.DensityPeakClustering(sRes, dRes, hCfg);
end