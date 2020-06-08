function loadFiles(obj)
    %LOADFILES Load results struct
    if obj.isError
        error(obj.errMsg);
    end

    if ~exist(obj.hCfg.resFile, 'file')
        error('%s does not exist', obj.hCfg.resFile);
        return;
    end

    obj.hCfg.updateLog('loadRes', sprintf('Loading %s', obj.hCfg.resFile), 1, 0);
    try
        res_ = load(obj.hCfg.resFile);
    catch ME
        %% try to some hacks to get around windows path name length limitations, if that is the problem
        if ispc
            userDir = char(java.lang.System.getProperty('user.home'));
            filename_temp = [userDir,filesep,'tmp_jrclust.mat']; 
            failed = system(sprintf('copy "%s" "%s"',obj.hCfg.resFile,filename_temp));        
            if ~failed
                try
                    res_ = load(filename_temp);
                    delete(filename_temp);
                catch
                    warning('Failed to load %s: %s', ME.message);
                    return;
                end
            else
                warning('Failed to load %s: %s', ME.message);
                return;            
            end  
        else
            warning('Failed to load %s: %s', ME.message);
            return;
        end
    end
    obj.hCfg.updateLog('loadRes', sprintf('Finished loading %s', obj.hCfg.resFile), 0, 1);


    if isfield(res_, 'spikeTimes')
        % load spikesRaw
        if isfield(res_, 'rawShape') && ~isempty(res_.rawShape)
            obj.hCfg.updateLog('loadRaw', sprintf('Loading %s', obj.hCfg.rawFile), 1, 0);
            spikesRaw = readBin(obj.hCfg.rawFile, res_.rawShape, '*int16');
            obj.hCfg.updateLog('loadRaw', sprintf('Finished loading %s', obj.hCfg.rawFile), 0, 1);
        else
            spikesRaw = [];
        end

        % load spikesFilt
        if isfield(res_, 'filtShape') && ~isempty(res_.filtShape)
            obj.hCfg.updateLog('loadFilt', sprintf('Loading %s', obj.hCfg.filtFile), 1, 0);
            spikesFilt = readBin(obj.hCfg.filtFile, res_.filtShape, '*int16');
            obj.hCfg.updateLog('loadFilt', sprintf('Finished loading %s', obj.hCfg.filtFile), 0, 1);
        else
            spikesFilt = [];
        end

        % load spikeFeatures
        if isfield(res_, 'featuresShape') && ~isempty(res_.featuresShape)
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

        if isfield(res_, 'history') && iscell(res_.history) % old-style history
            res_.history = convertHistory(res_.history, res_.initialClustering, obj.hCfg);
        end

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
            hClust = jrclust.sort.TemplateClustering(obj.hCfg);
            fieldNames = fieldnames(res_);

            md = ?jrclust.sort.TemplateClustering;
            pl = md.PropertyList;
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                propMetadata = pl(strcmp(fn, {pl.Name}));
                if isprop(hClust, fn) && ~((propMetadata.Dependent && isempty(propMetadata.SetMethod)))
                    hClust.(fn) = res_.(fn);
                end
            end

            res_.hClust = hClust;
        elseif isfield(res_, 'spikeClusters')
            hClust = jrclust.sort.DensityPeakClustering(obj.hCfg);
            fieldNames = fieldnames(res_);

            md = ?jrclust.sort.DensityPeakClustering;
            pl = md.PropertyList;
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                propMetadata = pl(strcmp(fn, {pl.Name}));
                if isprop(hClust, fn) && ~((propMetadata.Dependent && isempty(propMetadata.SetMethod)))
                    hClust.(fn) = res_.(fn);
                end
            end

            res_.hClust = hClust;
        end

        if isfield(res_, 'hClust')
            if ~isempty(res_.hClust.inconsistentFields()) && obj.hCfg.getOr('autoRecover', 0)
                flag = res_.hClust.recover(1); % recover inconsistent data if needed
                successAppend = 'You should look through your data and ensure everything is correct, then save it.';
                failureAppend = 'You will probably experience problems curating your data.';
                msg = '';
                switch flag
                    case 2
                        msg = sprintf('Non-contiguous spike table found and corrected. %s', successAppend);
                        
                    case 1
                        msg = sprintf('Inconsistent fields found and corrected. %s', successAppend);

                    case 0
                        msg = sprintf('Clustering data in an inconsistent state and automatic recovery failed. Please post an issue on the GitHub issue tracker. %s', failureAppend);
                        
                    case -1
                        msg = sprintf('Automatic recovery canceled by the user but the clustering data is still in an inconsistent state. %s', failureAppend);
                end

                if ~isempty(msg)
                    jrclust.utils.qMsgBox(msg, 1, 1);
                end
            end

            res_.hClust.syncHistFile();
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
    %LOADBIN Load traces/features from binary file
    if exist(filename, 'file')
        fid = fopen(filename, 'r');
        binData = fread(fid, Inf, dataType);
        fclose(fid);
        binData = reshape(binData, binShape);
    else
        binData = [];
    end
end

function hClust = convertToNew(res, hCfg)
    %CONVERTHISTORY Convert old-style hClust to new-style
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

    hClust = jrclust.sort.DensityPeakClustering(hCfg, sRes, dRes);
end

function history = convertHistory(oldHistory, spikeClusters, hCfg)
    %CONVERTHISTORY Convert old-style (cell) history to the history file
    history = containers.Map('KeyType', 'int32', 'ValueType', 'char');
    history(1) = 'initial commit';

    fidHist = fopen(hCfg.histFile, 'w');
    fwrite(fidHist, int32(1), 'int32');
    fwrite(fidHist, int32(spikeClusters), 'int32');

    for j = 2:size(oldHistory, 1)
        op = oldHistory(j, :);

        switch op{2}
            case 'delete'
                deleted = op{3};
                nClustersOld = max(spikeClusters);
                spikeClusters(ismember(spikeClusters, deleted)) = 0;
                if ~isempty(obj.clusterCenters)
                    obj.clusterCenters(deleted) = [];
                end

                % shift clusters larger than deleted down by 1
                if numel(deleted) == 1
                    gtMask = (spikeClusters > deleted);
                    spikeClusters(gtMask) = spikeClusters(gtMask) - 1;
                else
                    keepMe = setdiff(1:nClustersOld, deleted);
                    nClustersNew = numel(keepMe);

                    good = (spikeClusters > 0);
                    mapFrom = zeros(1, nClustersOld);
                    mapFrom(keepMe) = 1:nClustersNew;

                    spikeClusters(good) = mapFrom(spikeClusters(good));
                end

            case 'merge'
                iCluster = op{3};
                jCluster = op{4};
                spikeClusters(spikeClusters == jCluster) = iCluster;

                % shift clusters larger than jCluster down by 1
                gtMask = (spikeClusters > jCluster);
                spikeClusters(gtMask) = spikeClusters(gtMask) - 1;

            case 'split'
                iCluster = op{3};
                retained = op{4};

                % shift clusters larger than iCluster up by 1 (make
                % room for splitted off cluster)
                gtMask = (spikeClusters > iCluster);
                spikeClusters(gtMask) = spikeClusters(gtMask) + 1;

                % take splitted off spikes and make a new cluster
                % of them
                iMask = (spikeClusters == iCluster);
                iSpikes = find(iMask);
                jSpikes = iSpikes(~ismember(iSpikes, retained)); % spikes to split off
                spikeClusters(jSpikes) = iCluster + 1;

            case 'partition'
                iCluster = op{3};
                assignPart = op{4};
                initialSpikes = find(spikeClusters == iCluster);

                for kSplit = 1:numel(assignPart)-1
                    if kSplit > 1
                        iSpikes = find(spikeClusters == iCluster);
                    else
                        iSpikes = initialSpikes;
                    end

                    splitOff = initialSpikes(assignPart{kSplit});
                    retain = iSpikes(~ismember(iSpikes, splitOff));

                    iSite = mode(obj.spikeSites(retain));
                    jSite = mode(obj.spikeSites(splitOff));

                    % swap retain and split-off spikes; always make splitOff the next
                    % cluster
                    if iSite > jSite
                        splitOff = retain;
                    end

                    % make room for new cluster
                    mask = (spikeClusters > iCluster + kSplit - 1);
                    spikeClusters(mask) = spikeClusters(mask) + 1;
                    spikeClusters(splitOff) = iCluster + kSplit;
                end

            otherwise
                diffs = oldHistory{j, 3};
                if ~isempty(diffs)
                    iDiffs = diffs(1, :);
                    sDiffs = diffs(2, :);
                    spikeClusters(iDiffs) = sDiffs;
                end

        end

        fwrite(fidHist, int32(j), 'int32');
        fwrite(fidHist, spikeClusters, 'int32');
        history(j) = op{2};
    end

    fclose(fidHist);
end