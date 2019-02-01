classdef DensityPeakClustering < jrclust.interfaces.Clustering
    %DENSITYPEAKCLUSTERING A Rodriguez-Laio clustering of spike data
    properties (Hidden, SetObservable)
        hCfg;               % Config object
    end

    %% DETECTION/CLUSTERING RESULTS
    properties (Hidden, SetAccess=protected, SetObservable)
        sRes;               % sorting results
        dRes;               % detection results
    end

    %% OLD-STYLE PROPERTIES, publicly gettable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        cviSpk_clu;         % => spikesByCluster
        csNote_clu;         % => clusterNotes
        delta;              % => spikeDelta
        icl;                % => clusterCenters
        mrWavCor;           % => simScore
        nClu;               % => nClusters
        nneigh;             % => spikeNeigh
        ordrho;             % => ordRho
        rho;                % => spikeRho
        tmrWav_clu;         % => meanWfGlobal
        tmrWav_spk_clu;     % => meanWfGlobal
        tmrWav_raw_clu;     % => meanWfGlobalRaw
        trWav_spk_clu;      % => meanWfLocal
        trWav_raw_clu;      % => meanWfLocalRaw
        tmrWav_raw_hi_clu;  % => meanWfRawHigh
        tmrWav_raw_lo_clu;  % => meanWfRawLow
        viClu;              % => spikeClusters
        viClu_auto;         % => initialClustering
        viSite_clu;         % => clusterSites
        viSite_min_clu;     % => unitPeakSites
        vnSite_clu;         % => nSitesOverThresh
        vnSpk_clu;          % => unitCount
        vrIsiRatio_clu;     % => unitISIRatio
        vrIsoDist_clu;      % => unitIsoDist
        vrLRatio_clu;       % => unitLRatio
        vrPosX_clu;         % => clusterCentroids(:, 1)
        vrPosY_clu;         % => clusterCentroids(:, 2)
        vrSnr_clu;          % => unitSNR
        vrVmin_clu;         % => unitPeaks
        vrVmin_uv_clu;      % => unitPeaksRaw
        vrVpp_clu;          % => unitVpp
        vrVpp_uv_clu;       % => unitVppRaw
        vrVrms_site;        % => siteRMS
    end

    %% NEW-STYLE CLUSTERING PROPERTIES
    properties (SetAccess=private, SetObservable)
        clusterCenters;     % cluster centers
        clusterCentroids;   % centroids of clusters on the probe
        unitCount;          % number of spikes per cluster
        clusterNotes;       % notes on clusters
        clusterSites;       % site on which spikes in this cluster most often occur
        editPos;            % current position in edit history
        history;            % cell array, log of merge/split/delete operations
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
        meanWfLocal;        % mean filtered waveforms for each cluster
        meanWfGlobal;       % mean filtered waveforms for each cluster over all sites
        meanWfLocalRaw;     % mean raw waveforms for each cluster
        meanWfGlobalRaw;    % mean raw waveforms for each cluster over all sites
        meanWfRawLow;       % mean raw waveforms for each cluster over all sites at a low point on the probe (for drift correction)
        meanWfRawHigh;      % mean raw waveforms for each cluster over all sites at a high point on the probe (for drift correction)
        simScore;           % waveform-based cluster similarity scores
    end

    % quality metrics
    properties (SetAccess=private, SetObservable)
        nSitesOverThresh;   % number of sites exceeding the detection threshold, per cluster
        siteRMS;            % site-wise 
        unitPeaks;          % minimum voltage of mean filtered waveforms at peak site, per cluster
        unitPeaksRaw;       % minimum voltage (uV) of mean raw waveforms at peak site, per cluster
        unitPeakSites;      % sites on which unitPeaks occur
        unitVpp;            % peak-to-peak voltage of filtered waveforms at peak site, per cluster
        unitVppRaw;         % peak-to-peak voltage of raw waveforms at peak site, per cluster
        unitISIRatio;       % inter-spike interval ratio #(ISI <= 2ms)/#(ISI <= 20ms), per cluster
        unitIsoDist;        % isolation distance
        unitLRatio;         % L-ratio
        unitSNR;            % signal-to-noise ratio at peak site (peak/RMS)
    end

    % computed from other values, but only on set
    properties (SetAccess=private, Transient)
        nClusters;          % number of clusters
    end

    % computed from other values
    properties (Dependent, Transient)
        nEdits;             % number of edits made to initial clustering
    end

    %% SORTING RESULTS, IMMUTABLE
    properties (Dependent, Transient)
        initialClustering;  % initial assignment of spikes to cluster
        ordRho;             % spike-wise index ordered by density (DPCLUS)
        rhoCuts;            % site-wise distance cutoff values for computing density (DPCLUS)
        spikeDelta;         % spike-wise distance to nearest neighbor of higher density (DPCLUS)
        spikeNeigh;         % nearest neighbor of higher density (DPCLUS)
        spikeRho;           % spike-wise density (DPCLUS)
    end

    %% DETECTION RESULTS
    properties (Dependent, Transient)
        siteThresh;         % sitewise detection threshold
        spikeAmps;          % amplitudes of detected spikes
        spikeFeatures;      % features which were clustered
        spikePositions;     % positions on the probe at which spikes are detected
        spikesBySite;       % aggregate of spike indices by site
        spikesBySite2;      % aggregate of secondary spike indices by site
        spikesBySite3;      % aggregate of tertiary spike indices by site
        spikesFilt;         % filtered spike windows
        spikeSites;         % sites on which spikes occur
        spikeSites2;        % secondary sites on which spikes occur
        spikesRaw;          % raw spike windows
        spikeTimes;         % times at which spikes occurred
    end

    %% LIFECYCLE
    methods
        function obj = DensityPeakClustering(sRes, dRes, hCfg)
            obj.dRes = dRes;
            obj.hCfg = hCfg;
            obj.history = cell(0, 4);
            isImport = obj.tryImport(sRes);

            if ~isImport
                obj.sRes = sRes;
                obj.spikeClusters = obj.initialClustering;
            end

            % these fields are mutable so we need to store copies in obj
            if isfield(sRes, 'clusterCenters')
                obj.clusterCenters = sRes.clusterCenters;
            else
                obj.clusterCenters = [];
            end
            if isfield(sRes, 'clusterCentroids')
                obj.clusterCentroids = sRes.clusterCentroids;
            else
                obj.clusterCentroids = [];
            end

            if ~isImport
                obj.clearNotes();
                obj.refresh(1, []);
                commitMsg = sprintf('%s;initial commit', datestr(now, 31));
                obj.commit(commitMsg);
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function selfSim = computeSelfSim(obj, iCluster)
            %COMPUTESELFSIM Get similarity between bottom and top half (vpp-wise) of cluster
            if nargin < 2
                iCluster = [];
            end

            if isempty(iCluster)
                if obj.hCfg.verbose
                    fprintf('Computing self correlation\n\t');
                    t1 = tic;
                end

                selfSim = zeros(1, obj.nClusters);
                for iCluster = 1:obj.nClusters
                    selfSim(iCluster) = doComputeSelfSim(obj, iCluster);
                    if obj.hCfg.verbose
                        fprintf('.');
                    end
                end

                if obj.hCfg.verbose
                    fprintf('\n\ttook %0.1fs\n', toc(t1));
                end
            else
                selfSim = doComputeSelfSim(obj, iCluster);
            end
        end

        function [sites1, sites2, sites3] = getSecondaryPeaks(obj)
            %GETSECONDARYPEAKS
            minVals = squeeze(min(obj.meanWfGlobal) - obj.meanWfGlobal(1, :, :));

            [~, sites1] = min(minVals);

            minVals(sub2ind(size(minVals), sites1, 1:numel(sites1))) = 0;
            [~, sites2] = min(minVals);

            minVals(sub2ind(size(minVals), sites2, 1:numel(sites2))) = 0;
            [~, sites3] = min(minVals);
        end

        function nMerged = mergeBySim(obj)
            %MERGEBYSIM Automatically merge clusters by similarity score
            simScore_ = obj.simScore;
            nClusters_ = size(simScore_, 2);

            % identify clusters to remove, update and same (no change), disjoint sets
            simScore_(tril(true(nClusters_))) = 0; % ignore bottom half
            [scoresMax, mapTo] = max(simScore_);

            % keep clusters whose maximum similarity to some other cluster
            % is LESS than our threshold
            keepMe_ = scoresMax < obj.hCfg.maxUnitSim;

            if all(keepMe_)
                nMerged = 0;
                return;
            end

            minScore = min(scoresMax(~keepMe_));
            keepMe = find(keepMe_);
            mapTo(keepMe_) = keepMe; % map units to keep to themselves
 
            keepMe = setdiff(keepMe, mapTo(~keepMe_));
            removeMe = setdiff(1:nClusters_, mapTo);
            updateMe = setdiff(setdiff(1:nClusters_, keepMe), removeMe);

            spikeClusters_ = obj.spikeClusters;
            good = spikeClusters_ > 0;
            spikeClusters_(good) = int32(mapTo(spikeClusters_(good))); % translate cluster number

            obj.subsetFields(union(keepMe, updateMe));
            [~, ~, obj.spikeClusters(good)] = unique(spikeClusters_(good)); % remap to 1:nNewClusters

            updateMe = arrayfun(@(i) find(union(keepMe, updateMe) == i), updateMe);
            obj.refresh(0, updateMe); % recount

            arrayfun(@obj.rmRefracSpikes, updateMe); % remove refrac spikes
            obj.removeEmptyClusters();

            % update cluster waveforms and distance
            obj.computeMeanWaveforms(updateMe);
            obj.computeWaveformSim(updateMe);

            nMerged = nClusters_ - obj.nClusters;
            if obj.hCfg.verbose
                fprintf('\tnClusters: %d->%d (%d merged, min score: %0.4f)\n', nClusters_, obj.nClusters, nMerged, minScore);
            end
        end

        function postOp(obj, updateMe)
            %POSTOP Call this after a merge or split operation
            % update counts, center sites, remove empty clusters
            obj.refresh(1, updateMe);

            if ~isempty(updateMe)
                arrayfun(@obj.rmRefracSpikes, updateMe);
                %obj.rmRefracSpikes(updateMe);

                % update cluster waveforms and distance
                obj.updateWaveforms(updateMe);

                % update cluster positions
                obj.computeCentroids(updateMe);

                % compute quality scores for new clusters
                obj.computeQualityScores(updateMe);
            end
        end

        function removeEmptyClusters(obj)
            %REMOVEEMPTYCLUSTERS Find and remove empty clusters
            keepClusters = obj.unitCount > 0;
            if all(keepClusters)
                return;
            end

            % subset all fields indexed by cluster
            obj.subsetFields(keepClusters);

            if min(obj.spikeClusters) < 1 % noise or garbage clusters
                obj.spikeClusters(obj.spikeClusters < 1) = 0;

                % renumber clusters 1:numel(unique(spikeClusters))
                [~, ~, obj.spikeClusters] = unique(obj.spikeClusters + 1);
                obj.spikeClusters = obj.spikeClusters - 1;
            else
                % renumber clusters 1:numel(unique(spikeClusters))
                [~, ~, obj.spikeClusters] = unique(obj.spikeClusters);
            end

            obj.spikeClusters = int32(obj.spikeClusters);
        end

        function subsetFields(obj, keepMe)
            %SELECT Subset all data fields, taking only those indices we
            %want to keep, prior to rearranging or deleting
            fieldNames = fieldnames(obj);
            vecFields = fieldNames(cellfun(@(fn) isvector(obj.(fn)) && numel(obj.(fn)) == obj.nClusters, fieldNames));
            for i = 1:numel(vecFields)
                fn = vecFields{i};
                fd = obj.(fn);
                obj.(fn) = fd(keepMe);
            end

            % subset matrix fields
            if ~isempty(obj.clusterCentroids) % nClusters x 2
                obj.clusterCentroids = obj.clusterCentroids(keepMe, :);
            end
            if ~isempty(obj.simScore) % nClusters x nClusters
                obj.simScore = obj.simScore(keepMe, keepMe);
            end

            % subset tensor fields
            tenFields = fieldNames(cellfun(@(fn) ndims(obj.(fn)) == 3 && size(obj.(fn), 3) == obj.nClusters, fieldNames));
            for i = 1:numel(tenFields)
                fn = tenFields{i};
                fd = obj.(fn);
                obj.(fn) = fd(:, :, keepMe);
            end
        end

        function success = tryImport(obj, sRes)
            reqFields = {'spikeClusters', 'spikeRho', 'spikeDelta', ...
                         'spikeNeigh', 'ordRho', 'clusterCenters', ...
                         'initialClustering', 'clusterNotes'};
            if all(ismember(reqFields, fieldnames(sRes)))
                obj.spikeClusters = double(sRes.spikeClusters);
                % object takes initalClustering from sRes.spikeClusters
                if isfield(sRes, 'initialClustering')
                    sRes.spikeClusters = double(sRes.initialClustering);
                end

                if ~isfield(sRes, 'ordRho')
                    [~, sRes.ordRho] = sort(obj.spikeRho, 'descend');
                end

                % mean waveforms
                meanFields = {'meanWfGlobal', 'meanWfGlobalRaw', 'meanWfLocal', ...
                              'meanWfLocalRaw', 'meanWfRawHigh', 'meanWfRawLow'};
                for i = 1:numel(meanFields)
                    mf = meanFields{i};
                    if isfield(sRes, mf)
                        obj.(mf) = sRes.(mf);
                    end
                end

                % cluster-wise summaries
                sumFields = {'unitCount', 'clusterSites', 'spikesByCluster'};
                for i = 1:numel(sumFields)
                    mf = sumFields{i};
                    if isfield(sRes, mf)
                        obj.(mf) = sRes.(mf);
                    end
                end

                % quality scores
                qualFields = {'nSitesOverThresh', 'siteRMS', 'unitISIRatio', ...
                              'unitIsoDist', 'unitLRatio', 'unitPeaks', ...
                              'unitPeaksRaw', 'unitSNR', 'unitVpp', 'unitVppRaw'};
                for i = 1:numel(qualFields)
                    mf = qualFields{i};
                    if isfield(sRes, mf)
                        obj.(mf) = sRes.(mf);
                    end
                end

                % sim score
                if isfield(sRes, 'simScore')
                    obj.simScore = sRes.simScore;
                else
                    obj.simScore = doComputeWaveformSim(obj, []);
                end

                if isfield(sRes, 'clusterNotes')
                    obj.clusterNotes = sRes.clusterNotes;
                else
                    obj.clearNotes();
                end

                obj.sRes = sRes;
                % any mean waveform fields empty? recompute
                if any(cellfun(@(f) isempty(obj.(f)), meanFields))
                    obj.updateWaveforms([]);
                end

                % any summary fields empty? refresh
                if any(cellfun(@(f) isempty(obj.(f)), sumFields))
                    obj.refresh(1, []);
                end

                % any quality scores empty? recompute
                if any(cellfun(@(f) isempty(obj.(f)), qualFields))
                    obj.computeQualityScores();
                end
                commitMsg = sprintf('%s;initial import', datestr(now, 31));
                obj.commit(commitMsg);

                success = 1;
            else
                success = 0;
            end
        end
    end

    %% USER METHODS
    methods
        function addNote(obj, iCluster, note)
            %ADDNOTE Annotate a cluster
            if iCluster < 1 || iCluster > obj.nClusters
                return;
            end

            obj.clusterNotes{iCluster} = note;
        end

        function success = autoMerge(obj, maxUnitSim)
            %AUTOMERGE Automatically merge clusters
            success = 0;

            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            if nargin > 1
                try
                    obj.hCfg.setTemporaryParams('maxUnitSim', maxUnitSim);
                catch ME
                    warning('autoMerge aborted: %s', ME.message);
                    return;
                end
            end

            obj.refresh(1, []);
            obj.orderClusters('clusterSites');
            obj.clearNotes();

            obj.rmOutlierSpikes();
            obj.doWaveformMerge();
            obj.refresh(1, []);
            obj.orderClusters('clusterSites');
            obj.updateWaveforms();

            obj.computeCentroids();
            obj.clearNotes();
            obj.computeQualityScores();
            commitMsg = sprintf('%s;autoMerge', datestr(now, 31));
            obj.commit(commitMsg);

            if nargin == 2
                obj.hCfg.resetTemporaryParams('maxUnitSim');
            end

            success = 1;
        end

        function clearNotes(obj)
            %CLEARNOTES Remove all cluster notes
            obj.clusterNotes = arrayfun(@(~) '', 1:obj.nClusters, 'UniformOutput', 0);
        end

        function rates = firingRates(obj, clusters, nSamples)
            %FIRINGRATES Compute firing rates for specified clusters
            if nargin < 2
                clusters = 1:obj.nClusters;
            end
            if nargin < 3 || isempty(nSamples)
                nSamples = round(obj.hCfg.recDurationSec()*obj.hCfg.frSampleRate);
            end

            nFilt = round(obj.hCfg.frSampleRate * obj.hCfg.frPeriod / 2);
            if strcmp(obj.hCfg.frFilterShape, 'triangle')
                filtKernel = ([1:nFilt, nFilt-1:-1:1]'/nFilt*2/obj.hCfg.frPeriod);
            elseif strctmp(obj.hCfg.frFilterShape, 'rectangle')
                filtKernel = (ones(nFilt*2, 1) / obj.hCfg.frPeriod);
            end % switch

            filtKernel = single(filtKernel);

            rates = zeros([nSamples, numel(clusters)], 'single');
            for iiCluster = 1:numel(clusters)
                iCluster = clusters(iiCluster);
                clusterSpikes = obj.spikesByCluster{iCluster};
                clusterTimes = obj.spikeTimes(clusterSpikes);

                dsTimes = round(obj.hCfg.frSampleRate*double(clusterTimes)/obj.hCfg.sampleRate);
                dsTimes = max(min(dsTimes, nSamples), 1);

                rates(dsTimes, iiCluster) = 1;
                rates(:, iiCluster) = conv(rates(:, iiCluster), filtKernel, 'same');
            end
        end

        function commit(obj, msg)
            %COMMIT Commit a modification of clustering to history log
            %   Message format: parts are separated by a semicolon, no spaces
            %   1: datestr(timestamp, 31)
            %   2: operation (initial commit/import; delete; merge; split)
            %   3: cluster(s) operated on:  for delete, all clusters deleted, comma-separated
            %                               for merge, lower-valued cluster
            %                               for split, cluster splitted
            %   4: operation-specific data: for delete, empty
            %                               for merge, higher-valued cluster
            %                               for split, spikes retained, comma-separated
            %   Example: 2018-12-26 13:02:17;merge;4;5 denotes cluster 5
            %   merged into cluster 4
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            % check for consistency before committing
            ic = obj.selfConsistent();
            if ~isempty(ic)
                wmsg = strjoin(ic, '\n\t');
                warning('Cluster data inconsistent after previous operations:\n\t%s', wmsg);
                obj.editSeek(obj.nEdits);
                return;
            end

            % split commit message by commas
            msgParts = strsplit(lower(msg), ';');
            if numel(msgParts) < 2
                warning('malformed commit message: %s', msg);
                return
            end

            obj.history(end+1, 1:2) = msgParts(1:2);
            if startsWith(msgParts{2}, 'delete')
                % multiple clusters can be deleted at once
                deleted = cellfun(@(x) str2double(x), strsplit(msgParts{3}, ','));
                obj.history{end, 3} = deleted;
            elseif startsWith(msgParts{2}, 'merge')
                iCluster = str2double(msgParts{3});
                jCluster = str2double(msgParts{4});
                obj.history{end, 3} = min(iCluster, jCluster);
                obj.history{end, 4} = max(iCluster, jCluster);
            elseif startsWith(msgParts{2}, 'split')
                obj.history{end, 3} = str2double(msgParts{3});
                retained = cellfun(@(x) str2double(x), strsplit(msgParts{4}, ','));
                obj.history{end, 4} = retained;
            elseif any(obj.spikeClusters ~= obj.initialClustering) % store diff from previous clustering
                spikeClusters_ = obj.initialClustering;

                for j = 1:size(obj.history, 1)-1
                    op = obj.history(j, :);
                    if strcmp(op{2}, 'delete')
                        deleted = op{3};
                        spikeClusters_(ismember(spikeClusters_, deleted)) = 0;
                    elseif strcmp(op{2}, 'merge')
                        iCluster = op{3};
                        jCluster = op{4};
                        spikeClusters_(spikeClusters_ == jCluster) = iCluster;

                        % shift clusters larger than jCluster down by 1
                        gtMask = (spikeClusters_ > jCluster);
                        spikeClusters_(gtMask) = spikeClusters_(gtMask) - 1;
                    elseif strcmp(op{2}, 'split')
                        iCluster = op{3};
                        retained = op{4};

                        % shift clusters larger than iCluster up by 1 (make
                        % room for splitted off cluster)
                        gtMask = (spikeClusters_ > iCluster);
                        spikeClusters_(gtMask) = spikeClusters_(gtMask) + 1;

                        % take splitted off spikes and make a new cluster
                        % of them
                        newMask = (spikeClusters_ == iCluster & ~ismember(spikeClusters_, retained));
                        spikeClusters_(newMask) = iCluster + 1;
                    else
                        diffs = obj.history{j, 3};
                        if isempty(diffs)
                            continue;
                        end

                        iDiffs = diffs(1, :);
                        sDiffs = diffs(2, :);
                        spikeClusters_(iDiffs) = sDiffs;
                    end
                end % now caught up to the current state, can get diff

                iDiffs = find(obj.spikeClusters ~= spikeClusters_);
                sDiffs = obj.spikeClusters(iDiffs);
                obj.history{end, 3} = [iDiffs'; sDiffs'];
            end

            obj.editPos = size(obj.history, 1) - 1; % 0 for initial commmit, etc.
        end

        function computeCentroids(obj, updateMe)
            %COMPUTEPOSITIONS determine cluster position from spike position
            if nargin < 2 || isempty(obj.clusterCentroids)
                updateMe = [];
            end

            if isempty(updateMe)
                centroids = zeros(obj.nClusters, 2);
                clusters_ = 1:obj.nClusters;
            else % selective update
                centroids = obj.clusterCentroids;
                clusters_ = updateMe(:)';
            end

            featureSites = 1:obj.hCfg.nSitesEvt;

            for iCluster = clusters_
                [clusterSpikes, neighbors] = subsampleCenteredSpikes(obj, iCluster);
                if isempty(clusterSpikes)
                    continue;
                end

                neighbors = neighbors(1:end-obj.hCfg.nSitesExcl);
                featureWeights = squeeze(obj.spikeFeatures(featureSites, 1, clusterSpikes));
                neighborLoc = single(obj.hCfg.siteLoc(neighbors, :)); % position on probe

                centroids(iCluster, 1) = median(getWeightedLoc(featureWeights, neighborLoc(:, 1)));
                centroids(iCluster, 2) = median(getWeightedLoc(featureWeights, neighborLoc(:, 2)));
            end

            obj.clusterCentroids = centroids;
        end

        function computeMeanWaveforms(obj, updateMe, useRaw)
            %COMPUTEMEANWAVEFORMS Compute the mean waveform for each cluster
            if nargin < 2
                updateMe = [];
            end
            if nargin < 3
                useRaw = 1;
            end

            results = doComputeMeanWaveforms(obj, updateMe, useRaw);
            
            % collect computed values
            obj.unitPeaks = results.unitPeaks;
            obj.unitPeakSites = results.unitPeakSites;
            obj.meanWfLocal = results.meanWfLocal;
            obj.meanWfGlobal = results.meanWfGlobal;
            obj.meanWfLocalRaw = results.meanWfLocalRaw;
            obj.meanWfGlobalRaw = results.meanWfGlobalRaw;
            obj.meanWfRawLow = results.meanWfRawLow;
            obj.meanWfRawHigh = results.meanWfRawHigh;

        end

        function computeWaveformSim(obj, updateMe)
            %WAVEFORMSIM Compute waveform-based similarity scores for all clusters
            if nargin < 2 || isempty(obj.simScore)
                updateMe = [];
            end

            if isempty(obj.meanWfGlobal)
                return;
            end

            obj.simScore = doComputeWaveformSim(obj, updateMe);
            obj.simScore = jrclust.utils.setDiag(obj.simScore, obj.computeSelfSim());
        end

        function computeQualityScores(obj, updateMe)
            %COMPUTEQUALITYSCORES Get cluster quality scores
            if nargin < 2
                updateMe = [];
            end
            scores = doComputeQualityScores(obj, updateMe);
            obj.nSitesOverThresh = scores.nSitesOverThresh;
            obj.siteRMS = scores.siteRMS;
            obj.unitISIRatio = scores.unitISIRatio;
            obj.unitIsoDist = scores.unitIsoDist;
            obj.unitLRatio = scores.unitLRatio;
            obj.unitPeaksRaw = scores.unitPeaksRaw; % unitPeaks set elsewhere
            obj.unitSNR = scores.unitSNR;
            obj.unitVpp = scores.unitVpp;
            obj.unitVppRaw = scores.unitVppRaw;
        end

        function success = deleteClusters(obj, deleteMe)
            %DELETECLUSTERS Delete clusters
            success = 1;

            clustersBak = obj.spikeClusters; % in case we need to restore

            nClustersOld = numel(obj.spikesByCluster);
            keepMe = setdiff(1:nClustersOld, deleteMe);
            obj.subsetFields(keepMe);

            garbageCluster = min(obj.spikeClusters) - 1;
            if garbageCluster == 0 % reserved for noise
                garbageCluster = -1;
            end

            deleteSpikes = ismember(obj.spikeClusters, deleteMe);
            obj.spikeClusters(deleteSpikes) = garbageCluster;
            nClusters_ = numel(keepMe);

            good = (obj.spikeClusters > 0);
            mapFrom = zeros(1, nClustersOld);
            mapFrom(keepMe) = 1:nClusters_;

            obj.spikeClusters(good) = mapFrom(obj.spikeClusters(good));

            % recompute similarity scores
%             if isfield(obj, 'mrSim_clu')
%                 obj = sim_score_(obj);
%             end

            if ~isempty(obj.selfConsistent)
                warning('Cluster data is inconsistent after deleting %d', deleteMe);
                obj.spikeClusters = clustersBak;
                success = 0;
                obj.refresh(1, []);
            end
        end

        function doWaveformMerge(obj)
            %DOWAVEFORMMERGE Merge clusters by waveform-based similarity scores
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            obj.computeMeanWaveforms();
            obj.computeWaveformSim();

            for iRepeat = 1:obj.hCfg.nPassesMerge % single-pass vs dual-pass correction
                nMerged = obj.mergeBySim();
                if nMerged < 1
                    break;
                end
            end
        end

        function editSeek(obj, seekTo)
            %EDITSEEK Seek back and forth to position in history
            if seekTo < 0 || seekTo > obj.nEdits
                return;
            end

            % restore initial state and fast forward
            spikeClusters_ = obj.initialClustering;
            if isfield(obj.sRes, 'clusterCenters')
                obj.clusterCenters = obj.sRes.clusterCenters;
            else
                obj.clusterCenters = [];
            end

            for j = 1:seekTo+1
                op = obj.history(j, :);
                if strcmp(op{2}, 'delete')
                    deleted = op{3};
                    spikeClusters_(ismember(spikeClusters_, deleted)) = 0;
                elseif strcmp(op{2}, 'merge')
                    iCluster = op{3};
                    jCluster = op{4};
                    spikeClusters_(spikeClusters_ == jCluster) = iCluster;

                    % shift clusters larger than jCluster down by 1
                    gtMask = (spikeClusters_ > jCluster);
                    spikeClusters_(gtMask) = spikeClusters_(gtMask) - 1;
                elseif strcmp(op{2}, 'split')
                    iCluster = op{3};
                    retained = op{4};

                    % shift clusters larger than iCluster up by 1 (make
                    % room for splitted off cluster)
                    gtMask = (spikeClusters_ > iCluster);
                    spikeClusters_(gtMask) = spikeClusters_(gtMask) + 1;

                    % take splitted off spikes and make a new cluster
                    % of them
                    newMask = (spikeClusters_ == iCluster & ~ismember(spikeClusters_, retained));
                    spikeClusters_(newMask) = iCluster + 1;
                else
                    diffs = obj.history{j, 3};
                    if isempty(diffs)
                        continue;
                    end

                    iDiffs = diffs(1, :);
                    sDiffs = diffs(2, :);
                    spikeClusters_(iDiffs) = sDiffs;
                end
            end
            obj.spikeClusters = spikeClusters_;

            obj.refresh(0, []);
            if ~isempty(obj.meanWfGlobal) % only compute mean waveforms if we already had them
                obj.updateWaveforms();
            end

            if ~isempty(obj.unitVpp)
                obj.computeQualityScores();
            end
        end

        function success = exportToCSV(obj, zeroIndex)
            %EXPORTTOCSV Export spike times, sites, and clusters to CSV
            if nargin < 2
                zeroIndex = 0;
            end

            spikeTimes_ = double(obj.spikeTimes) / obj.hCfg.sampleRate;
            spikeSites_ = double(obj.spikeSites) - double(zeroIndex);

            filename = jrclust.utils.subsExt(obj.hCfg.configFile, '.csv');

            % write header
            try
                fid = fopen(filename, 'w');
                fprintf(fid, 'spikeTimes,spikeClusters,spikeSites\n');
                fclose(fid);
            catch ME
                warning('Failed to export: %s', ME.message);
                success = 0;
                return;
            end

            % write values
            dlmwrite(filename, [spikeTimes_(:), double(obj.spikeClusters(:)), spikeSites_(:)], 'precision', 9, '-append');

            if obj.hCfg.verbose
                fprintf('Wrote to %s. Columns:\n', filename);
                fprintf('\tColumn 1: Spike time (s)\n');
                fprintf('\tColumn 2: Unit# (positive #: valid units, 0: noise cluster, negative #: deleted clusters)\n');
                fprintf('\tColumn 3: Site# (starts with 1)\n');
            end

            success = 1;
        end

        function success = exportQualityScores(obj, zeroIndex, fGui)
            %EXPORTQUALITYSCORES Export cluster quality scores to CSV
            if nargin < 2
                zeroIndex = 0;
            end
            if nargin < 3
                fGui = 0;
            end
            
            ID = (1:obj.nClusters)';
            SNR = obj.unitSNR(:);
            centerSite = obj.clusterSites(:) - double(zeroIndex);
            nSpikes = obj.unitCount(:);
            xPos = obj.clusterCentroids(:, 1);
            yPos = obj.clusterCentroids(:, 2);
            uVmin = obj.unitPeaksRaw(:);
            uVpp = obj.unitVppRaw(:);
            IsoDist = obj.unitIsoDist(:);
            LRatio = obj.unitLRatio(:);
            ISIRatio = obj.unitISIRatio(:);
            note = obj.clusterNotes(:);

            filename = jrclust.utils.subsExt(obj.hCfg.configFile, '_quality.csv');

            try
                table_ = table(ID, SNR, centerSite, nSpikes, xPos, yPos, uVmin, uVpp, IsoDist, LRatio, ISIRatio, note);
                writetable(table_, filename);
            catch ME
                warning('Failed to export: %s', ME.message);
                success = 0;
                return;
            end

            if obj.hCfg.verbose
                disp(table_);
                helpText = {sprintf('Wrote to %s. Columns:', filename), ...
                            sprintf('\tColumn 1: ID: Unit ID'), ...
                            sprintf('\tColumn 2: SNR: |Vp/Vrms|; Vp: negative peak amplitude of the peak site; Vrms: SD of the Gaussian noise (estimated from MAD)'), ...
                            sprintf('\tColumn 3: centerSite: Peak site number which contains the most negative peak amplitude'), ...
                            sprintf('\tColumn 4: nSpikes: Number of spikes'), ...
                            sprintf('\tColumn 5: xPos: x position (width dimension) center-of-mass'), ...
                            sprintf('\tColumn 6: yPos: y position (depth dimension) center-of-mass, referenced from the tip'), ...
                            sprintf('\tColumn 7: uVmin: Min. voltage (uV) of the mean raw waveforms at the peak site (microvolts)'), ...
                            sprintf('\tColumn 8: uVpp: peak-to-peak voltage (microvolts)'), ...
                            sprintf('\tColumn 9: IsoDist: Isolation distance quality metric'), ...
                            sprintf('\tColumn 10: LRatio: L-ratio quality metric'), ...
                            sprintf('\tColumn 11: ISIRatio: ISI-ratio quality metric'), ...
                            sprintf('\tColumn 12: note: user comments')};

                cellfun(@(x) fprintf('%s\n', x), helpText);
                if fGui
                    jrclust.utils.qMsgBox(helpText);
                end
            end

            success = 1;
        end

        function uInfo = exportUnitInfo(obj, iCluster)
            %EXPORTUNITINFO Get all data pertinent to a cluster
            if ~ismember(iCluster, 1:obj.nClusters)
                uInfo = [];
                return;
            end

            iSite = obj.clusterSites(iCluster);
            iNeighbors = obj.hCfg.siteNeighbors(:, iSite);

            pos = sprintf('Unit %d (x,y):(%0.1f, %0.1f)[pix]', iCluster, obj.clusterCentroids/obj.hCfg.umPerPix);

            % subsample some (raw or filtered) waveforms
            iSubset = jrclust.utils.subsample(obj.getCenteredSpikes(iCluster), obj.hCfg.nSpikesFigWav);
            if obj.hCfg.showRaw
                meanWf = obj.meanWfGlobalRaw(:, iNeighbors, iCluster);
                sampleWf = jrclust.utils.rawTouV(obj.spikesRaw(:, :, iSubset), obj.hCfg);
                sampleWf = jrclust.filters.fftLowpass(sampleWf, obj.hCfg.getOr('fc_spkwav_show', []), obj.hCfg.sampleRate);
            else
                meanWf = obj.meanWfGlobal(:,iNeighbors,iCluster);
                sampleWf = jrclust.utils.filtTouV(obj.spikesFilt(:, :, iSubset), obj.hCfg);
            end

            uInfo = struct('cluster', iCluster, ...
                           'xyPos', obj.clusterCentroids(iCluster), ...
                           'meanWf', meanWf, ...
                           'neighbors', iNeighbors, ...
                           'position', pos, ...
                           'sampleWf', sampleWf);
            if ~isempty(obj.unitLRatio)
                uInfo.LRatio = obj.unitLRatio(iCluster);
            end
            if ~isempty(obj.unitISIRatio)
                uInfo.ISIRatio = obj.unitISIRatio(iCluster);
            end
            if ~isempty(obj.unitIsoDist)
                uInfo.IsoDist = obj.unitIsoDist(iCluster);
            end
            if ~isempty(obj.unitSNR)
                uInfo.SNR = obj.unitSNR(iCluster);
            end
            if ~isempty(obj.unitPeaksRaw)
                uInfo.peaksRaw = obj.unitPeaksRaw(iCluster);
            end
            if ~isempty(obj.unitVppRaw)
                uInfo.vpp = obj.unitVppRaw(iCluster);
            end
        end

        function [centeredSpikes, whichCentered] = getCenteredSpikes(obj, iCluster)
            %GETCENTEREDSPIKES Return subset of spikes which occur on center site of cluster
            iClusterSite = obj.clusterSites(iCluster);
            centeredSpikes = obj.spikesByCluster{iCluster};
            clusterSites_ = obj.spikeSites(centeredSpikes);
            whichCentered = clusterSites_ == iClusterSite;
            centeredSpikes = centeredSpikes(whichCentered);
        end

        function success = mergeClusterPair(obj, iCluster, jCluster)
            %MERGECLUSTERPAIR Merge a pair of clusters
            success = 0;
            if iCluster == jCluster || numel(intersect([iCluster, jCluster], obj.spikeClusters)) < 2
                return;
            end

            % keep the smaller of the two and shift others left by one
            [iCluster, jCluster] = deal(min(iCluster, jCluster), max(iCluster, jCluster));

            clustersBak = obj.spikeClusters;
            obj.spikeClusters(obj.spikeClusters == jCluster) = iCluster;
            obj.nClusters = obj.nClusters + 1; % this is a hack

            % take the denser of the two centers, using delta as tiebreaker
            iCenter = obj.clusterCenters(iCluster);
            jCenter = obj.clusterCenters(jCluster);
            if obj.spikeRho(iCenter) < obj.spikeRho(jCenter)
                obj.clusterCenters(iCluster) = jCenter;
            elseif obj.spikeRho(iCenter) == obj.spikeRho(jCenter) && obj.spikeDelta(iCenter) < obj.spikeDelta(jCenter)
                obj.clusterCenters(iCluster) = jCenter;
            end % otherwise, keep iCenter

            if ~obj.deleteClusters(jCluster) % subsets fields
                warning('Failed to delete cluster %d', jCluster);
                obj.spikeClusters = clustersBak;
                return;
            end

            if ~isempty(obj.selfConsistent)
                warning('Cluster data is inconsistent after merging %d and %d', iCluster, jCluster);
                obj.spikeClusters = clustersBak;
                success = 0;
                obj.refresh(1, []);
            else
                success = 1;
                obj.postOp(iCluster);
                obj.orderClusters('clusterSites');
            end
        end

        function inconsistencies = selfConsistent(obj)
            %SELFCONSISTENT Check all fields have the correct sizes
            inconsistencies = {};
            % vector fields
            if ~isempty(obj.clusterCenters) && numel(obj.clusterCenters) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterCenters: expected %d, actual %d', obj.nClusters, numel(obj.clusterCenters));
            end
            if ~isempty(obj.unitCount) && numel(obj.unitCount) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitCount: expected %d, actual %d', obj.nClusters, numel(obj.unitCount));
            end
            if ~isempty(obj.clusterSites) && numel(obj.clusterSites) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterSites: expected %d, actual %d', obj.nClusters, numel(obj.clusterSites));
            end
            if ~isempty(obj.nSitesOverThresh) && numel(obj.nSitesOverThresh) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('nSitesOverThresh: expected %d, actual %d', obj.nClusters, numel(obj.nSitesOverThresh));
            end
            if ~isempty(obj.unitISIRatio) && numel(obj.unitISIRatio) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitISIRatio: expected %d, actual %d', obj.nClusters, numel(obj.unitISIRatio));
            end
            if ~isempty(obj.unitIsoDist) && numel(obj.unitIsoDist) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitIsoDist: expected %d, actual %d', obj.nClusters, numel(obj.unitIsoDist));
            end
            if ~isempty(obj.unitLRatio) && numel(obj.unitLRatio) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitLRatio: expected %d, actual %d', obj.nClusters, numel(obj.unitLRatio));
            end
            if ~isempty(obj.unitPeakSites) && numel(obj.unitPeakSites) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitPeakSites: expected %d, actual %d', obj.nClusters, numel(obj.unitPeakSites));
            end
            if ~isempty(obj.unitPeaks) && numel(obj.unitPeaks) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitPeaks: expected %d, actual %d', obj.nClusters, numel(obj.unitPeaks));
            end
            if ~isempty(obj.unitPeaksRaw) && numel(obj.unitPeaksRaw) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitPeaksRaw: expected %d, actual %d', obj.nClusters, numel(obj.unitPeaksRaw));
            end
            if ~isempty(obj.unitSNR) && numel(obj.unitSNR) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitSNR: expected %d, actual %d', obj.nClusters, numel(obj.unitSNR));
            end
            if ~isempty(obj.unitVpp) && numel(obj.unitVpp) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitVpp: expected %d, actual %d', obj.nClusters, numel(obj.unitVpp));
            end
            if ~isempty(obj.unitVppRaw) && numel(obj.unitVppRaw) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('unitVppRaw: expected %d, actual %d', obj.nClusters, numel(obj.unitVppRaw));
            end

            % matrix fields
            if ~isempty(obj.clusterCentroids) && size(obj.clusterCentroids, 1) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterCentroids: expected %d, actual %d', obj.nClusters, size(obj.clusterCentroids, 1));
            end
            if ~isempty(obj.simScore) && ~all(size(obj.simScore) == obj.nClusters)
                inconsistencies{end+1} = sprintf('simScore: expected %dx%d, actual %dx%d', obj.nClusters, obj.nClusters, size(obj.simScore, 1), size(obj.simScore, 2));
            end

            % tensor fields
            if ~isempty(obj.meanWfGlobal) && size(obj.meanWfGlobal, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfGlobal: expected %d, actual %d', obj.nClusters, size(obj.meanWfGlobal, 3));
            end
            if ~isempty(obj.meanWfGlobalRaw) && size(obj.meanWfGlobalRaw, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfGlobalRaw: expected %d, actual %d', obj.nClusters, size(obj.meanWfGlobalRaw, 3));
            end
            if ~isempty(obj.meanWfLocal) && size(obj.meanWfLocal, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfLocal: expected %d, actual %d', obj.nClusters, size(obj.meanWfLocal, 3));
            end
            if ~isempty(obj.meanWfLocalRaw) && size(obj.meanWfLocalRaw, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfLocalRaw: expected %d, actual %d', obj.nClusters, size(obj.meanWfLocalRaw, 3));
            end
            if ~isempty(obj.meanWfRawHigh) && size(obj.meanWfRawHigh, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfRawHigh: expected %d, actual %d', obj.nClusters, size(obj.meanWfRawHigh, 3));
            end
            if ~isempty(obj.meanWfRawLow) && size(obj.meanWfRawLow, 3) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('meanWfRawLow: expected %d, actual %d', obj.nClusters, size(obj.meanWfRawLow, 3));
            end

            % cell fields
            if ~isempty(obj.clusterNotes) && numel(obj.clusterNotes) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterNotes: expected %d, actual %d', obj.nClusters, numel(obj.clusterNotes));
            end
            if ~isempty(obj.spikesByCluster) && numel(obj.spikesByCluster) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('spikesByCluster: expected %d, actual %d', obj.nClusters, numel(obj.spikesByCluster));
            end

            % clusterCenters
            if ~isempty(obj.clusterCenters) && numel(unique(obj.clusterCenters)) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterCenters: expected %d unique, actual %d', obj.nClusters, numel(unique(obj.clusterCenters)));
            end
        end

        function orderClusters(obj, by)
            %ORDERCLUSTERS Arrange cluster ID numbers by some criterion
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            if nargin < 2 || isempty(by) || ~isprop(obj, by)
                by = 'clusterSites';
            end

            if strcmpi(by, 'Y + X') && ~isempty(obj.clusterCentroids)
                [~, argsort] = sort(sum(obj.clusterCentroids, 2), 'ascend');
            elseif isprop(obj, by)
                [~, argsort] = sort(obj.(by), 'ascend');
            end

            %obj.spikeClusters = mapIndex_(obj.spikeClusters, argsort);
            map(argsort) = 1:numel(argsort);
            mask = (obj.spikeClusters > 0);
            obj.spikeClusters(mask) = map(obj.spikeClusters(mask)); % do not map zeros

            % reorder data fields
            obj.subsetFields(argsort);
%             obj.spikesByCluster = obj.spikesByCluster(argsort);
%             obj.clusterSites = obj.clusterSites(argsort);
%             obj.unitCount = obj.unitCount(argsort);
%             obj.clusterNotes = obj.clusterNotes(argsort);
%             if ~isempty(obj.clusterCentroids)
%                 obj.clusterCentroids = obj.clusterCentroids(argsort, :);
%             end
        end

        function reassign(obj, recompute)
            %REASSIGN Reassign clusters, e.g., after a change of parameters
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end
            if nargin < 2
                recompute = ~isempty(obj.meanWfGlobal); % don't recompute mean waveforms if we didn't have them already
            end

            obj.sRes = jrclust.cluster.densitypeaks.assignClusters(obj.dRes, obj.sRes, obj.hCfg);

            % these fields are mutable so we need to store copies in obj
            obj.spikeClusters = obj.initialClustering;
            if isfield(obj.sRes, 'clusterCenters')
                obj.clusterCenters = obj.sRes.clusterCenters;
            else
                obj.clusterCenters = [];
            end
            obj.clusterCentroids = [];

            obj.clearNotes();
            obj.refresh(1, []);

            if recompute
                obj.updateWaveforms();
                obj.computeCentroids();
                obj.computeQualityScores();
            end

            obj.orderClusters('clusterSites');
            obj.clearNotes();

            commitMsg = sprintf('%s;reassign', datestr(now, 31));
            obj.commit(commitMsg);
        end

        function refresh(obj, doRemoveEmpty, updateMe)
            %REFRESH Recount and store spikes by cluster, optionally removing empty clusters
            if nargin < 2
                doRemoveEmpty = 0;
            end
            if nargin < 3
                updateMe = [];
            end

            if isempty(updateMe)
                obj.spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), 1:obj.nClusters, 'UniformOutput', 0);
                obj.unitCount = cellfun(@numel, obj.spikesByCluster);
                if ~isempty(obj.spikeSites)
                    obj.clusterSites = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), 1:obj.nClusters));
                else
                    obj.clusterSites = [];
                end
            else
                obj.spikesByCluster(updateMe) = arrayfun(@(iC) find(obj.spikeClusters == iC), updateMe, 'UniformOutput', 0);
                obj.unitCount(updateMe) = cellfun(@numel, obj.spikesByCluster(updateMe));
                if ~isempty(obj.spikeSites)
                    obj.clusterSites(updateMe) = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), updateMe));
                else
                    obj.clusterSites = [];
                end
            end

            if doRemoveEmpty
                obj.removeEmptyClusters();
            end
        end

        function revert(obj, revertTo)
            %REVERT Delete history
            if revertTo < 0 || revertTo >= obj.nEdits
                return;
            end

            obj.editSeek(revertTo);
            obj.history(revertTo+2:end, :) = [];
            obj.removeEmptyClusters();
            obj.computeCentroids();
            obj.computeQualityScores();
            obj.clearNotes();
            obj.editPos = revertTo;
        end

        function rmOutlierSpikes(obj)
            %RMOUTLIERSPIKES Mahalanobis-distance based outlier removal
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            if obj.hCfg.outlierThresh == 0
                return;
            end

            % This is troubling:
            % Warning: Matrix is close to singular or badly scaled. Results may be inaccurate.
            % > In mahal (line 49)
            % TODO: investigate
            warning off;
            for iCluster = 1:obj.nClusters
                iSite = obj.clusterSites(iCluster);
                iSpikes = obj.spikesByCluster{iCluster};

                if isempty(iSpikes)
                    continue;
                end

                iFeatures = squeeze(obj.spikeFeatures(:, 1, iSpikes));

                % if iSite is secondary site for any spikes in this cluster, prefer
                % features computed from those over primary spike features
                if size(obj.spikeFeatures, 2) >= 2
                    iFeatures2 = squeeze(obj.spikeFeatures(:, 2, iSpikes));
                    onSite2 = find(obj.spikeSites2(iSpikes) == iSite);
                    iFeatures(:, onSite2) = iFeatures2(:, onSite2);
                end

                iFeatures = iFeatures';

                try % MAD transform of log self-Mahalanobis distance
                    lastwarn(''); % reset last warning to catch it

                    iDist = jrclust.utils.madScore(log(mahal(iFeatures, iFeatures)));
                    [wstr, wid] = lastwarn();
                    if strcmp(wid, 'MATLAB:nearlySingularMatrix')
                        error(wstr);
                    end
                catch
                    continue;
                end

                exclude = (iDist > obj.hCfg.outlierThresh);
                if any(exclude)
                    obj.spikesByCluster{iCluster} = iSpikes(~exclude);
                    obj.spikeClusters(iSpikes(exclude)) = 0; % classify as noise
                    obj.unitCount(iCluster) = numel(obj.spikesByCluster{iCluster});
                end

                fprintf('.');
            end
            warning on;
        end

        function nRemoved = rmRefracSpikes(obj, iCluster)
            %RMREFRACSPIKES remove refractory spikes
            if obj.nEdits ~= obj.editPos % not at tip of edit history, back out
                warning('cannot branch from history; use revert() first');
                return;
            end

            if nargin == 1 % recurse for each cluster
                nRemoved = 0;
                for iCluster_ = 1:obj.nClusters
                    nRemoved_ = obj.rmRefracSpikes(iCluster_);
                    nRemoved = nRemoved + nRemoved_;
                end

                return;
            else
                nRemoved = 0;
                nSkipRefrac = 4;

                try
                    clusterSpikes_ = obj.spikesByCluster{iCluster};
                catch
                    clusterSpikes_ = find(obj.spikeClusters == iCluster);
                end

                if isempty(clusterSpikes_)
                    return;
                end

                clusterTimes_ = obj.spikeTimes(clusterSpikes_);

                % removal loop
                keepMe = true(size(clusterTimes_));
                while (1)
                    iKeepMe = find(keepMe);

                    inRefrac = find(diff(clusterTimes_(keepMe)) < obj.hCfg.refracIntSamp) + 1;
                    if isempty(inRefrac)
                        break;
                    end

                    keepMe(iKeepMe(inRefrac(1:nSkipRefrac:end))) = 0;
                end

                nRemoved = sum(~keepMe);
                obj.spikeClusters(clusterSpikes_(~keepMe)) = 0; % assign to noise cluster

                obj.spikesByCluster{iCluster} = clusterSpikes_(keepMe);
                obj.unitCount(iCluster) = sum(keepMe);
            end
        end

        function [success, retained] = splitCluster(obj, iCluster, retained)
            %SPLITCLUSTER Split a cluster
            success = 0;
            if iCluster < 1 || iCluster > obj.nClusters
                return;
            end

            iSpikes = find(obj.spikeClusters == iCluster);
            if ~all(ismember(retained, iSpikes))
                return;
            end

            clustersBak = obj.spikeClusters;
            splitOff = iSpikes(~ismember(iSpikes, retained));

            % swap retained and splitOff if iSite > jSite
            iSite = mode(obj.spikeSites(retained));
            jSite = mode(obj.spikeSites(splitOff));
            if iSite > jSite
                [retained, splitOff] = deal(splitOff, retained);
            end

            fieldNames = fieldnames(obj);
            % augment vector fields
            vecFields = fieldNames(cellfun(@(fn) isvector(obj.(fn)) && numel(obj.(fn)) == obj.nClusters, fieldNames));
            for iField = 1:numel(vecFields)
                fn = vecFields{iField};
                fd = obj.(fn);

                if iscell(fd)
                    fd{end+1} = []; %#ok<*AGROW>
                else
                    fd(end+1) = 0;
                end
                obj.(fn) = fd;
            end

            % augment matrix fields
            if ~isempty(obj.clusterCentroids) % nClusters x 2
                obj.clusterCentroids(end+1, :) = zeros(1, 2, 'like', obj.clusterCentroids);
            end
            if ~isempty(obj.simScore) % nClusters x nClusters
                obj.simScore = [obj.simScore zeros(obj.nClusters, 1, 'like', obj.simScore); ...
                                zeros(1, obj.nClusters + 1, 'like', obj.simScore)];
            end

            % augment tensor fields
            tenFields = fieldNames(cellfun(@(fn) ndims(obj.(fn)) == 3 && size(obj.(fn), 3) == obj.nClusters, fieldNames));
            for iField = 1:numel(tenFields)
                fn = tenFields{iField};
                fd = obj.(fn);
                fd(:, :, end+1) = zeros(size(fd, 1), size(fd, 2), 'like', fd);
                obj.(fn) = fd;
            end

            % make room for new cluster
            mask = (obj.spikeClusters > iCluster);
            obj.spikeClusters(mask) = obj.spikeClusters(mask) + 1;
            obj.spikeClusters(splitOff) = iCluster + 1;

            % swap fields
            obj.subsetFields([1:iCluster obj.nClusters iCluster+1:obj.nClusters-1]);

            % get the maximally dense spike in new cluster and take it as a center
            [~, splitCenter] = max(obj.spikeRho(splitOff));
            obj.clusterCenters(iCluster + 1) = splitOff(splitCenter);
            % in case we took iCluster's center with us
            if obj.clusterCenters(iCluster) == obj.clusterCenters(iCluster + 1)
                [~, retCenter] = max(obj.spikeRho(retained));
                obj.clusterCenters(iCluster) = retained(retCenter);
            end

            if isempty(obj.selfConsistent())
                success = 1;
                obj.postOp([iCluster, iCluster + 1]);
                obj.orderClusters('clusterSites');
            else
                warning('Cluster data is inconsistent after splitting %d', iCluster);
                obj.spikeClusters = clustersBak;
                success = 0;
                obj.subsetFields([1:iCluster iCluster+1:obj.nClusters+1]);
                obj.refresh(1, []);
            end
        end

        function updateWaveforms(obj, updateMe)
            %UPDATEWAVEFORMS Update mean waveforms and sim scores
            if nargin < 2
                updateMe = [];
            end

            obj.computeMeanWaveforms(updateMe);
            obj.computeWaveformSim(updateMe);
        end
    end

    %% GETTERS/SETTERS
    methods
        function c = get.icl(obj)
            c = obj.clusterCenters;
        end

        % clusterNotes/csNote_clu
        function cn = get.csNote_clu(obj)
            cn = obj.clusterNotes;
        end

        % unitCount/vnSpk_clu
        function cc = get.vnSpk_clu(obj)
            cc = obj.unitCount;
        end

        % clusterSites/viSite_clu
        function cs = get.viSite_clu(obj)
            cs = obj.clusterSites;
        end

        % dRes
        function set.dRes(obj, dr)
            if isstruct(dr)
                assert(isfield(dr, 'spikeTimes'), 'spikeTimes must not be missing');
                assert(isfield(dr, 'spikeSites'), 'spikeSites must not be missing');
            elseif isobject(dr)
                assert(isprop(dr, 'spikeTimes'), 'spikeTimes must not be missing');
                assert(isprop(dr, 'spikeSites'), 'spikeSites must not be missing');
            else
                error('type not recognized: %s', class(dr));
            end

            obj.dRes = dr;
        end

        % hCfg
        function set.hCfg(obj, hc)
            failMsg = 'hCfg must be an object of type jrclust.Config';
            assert(isa(hc, 'jrclust.Config'), failMsg);
            obj.hCfg = hc;
        end

        % initialClustering/viClu_auto
        function ic = get.initialClustering(obj)
            if ~isempty(obj.sRes)
                ic = obj.sRes.spikeClusters;
            else
                ic = [];
            end
        end
        function ic = get.viClu_auto(obj)
            ic = obj.initialClustering;
        end

        % meanWfGlobal/tmrWav_spk_clu
        function mw = get.tmrWav_spk_clu(obj)
            mw = obj.meanWfGlobal;
        end

        % meanWfGlobalRaw/tmrWav_raw_clu
        function mw = get.tmrWav_raw_clu(obj)
            mw = obj.meanWfGlobalRaw;
        end

        % meanWfLocal/trWav_spk_clu
        function mw = get.trWav_spk_clu(obj)
            mw = obj.meanWfLocal;
        end

        % meanWfLocalRaw/trWav_raw_clu
        function mw = get.trWav_raw_clu(obj)
            mw = obj.meanWfLocalRaw;
        end

        % meanWfRawHigh/tmrWav_raw_hi_clu
        function mw = get.tmrWav_raw_hi_clu(obj)
            mw = obj.meanWfRawHigh;
        end

        % meanWfRawLow/tmrWav_raw_lo_clu
        function mw = get.tmrWav_raw_lo_clu(obj)
            mw = obj.meanWfRawLow;
        end

        % unitPeakSites/viSite_min_clu
        function ms = get.viSite_min_clu(obj)
            ms = obj.unitPeakSites;
        end

        % nClusters/nClu
        function nc = get.nClu(obj)
            nc = obj.nClusters;
        end

        % nEdits
        function ne = get.nEdits(obj)
            ne = size(obj.history, 1) - 1;
        end

        % ordRho/ordrho
        function or = get.ordRho(obj)
            if isfield(obj.sRes, 'ordRho')
                or = obj.sRes.ordRho;
            else
                or = [];
            end
        end
        function or = get.ordrho(obj)
            or = obj.ordRho;
        end

        % rhoCuts
        function rc = get.rhoCuts(obj)
            if isfield(obj.sRes, 'rhoCutSite')
                rc = obj.sRes.rhoCutSite;
            else
                rc = [];
            end
        end

        % simScore/mrWavCor
        function ss = get.mrWavCor(obj)
            ss = obj.simScore;
        end

        % nSitesOverThresh/vnSite_clu
        function so = get.vnSite_clu(obj)
            so = obj.nSitesOverThresh;
        end

        % siteThresh
        function st = get.siteThresh(obj)
            if isfield(obj.dRes, 'siteThresh') || isprop(obj.dRes, 'siteThresh')
                st = obj.dRes.siteThresh;
            else
                st = [];
            end
        end

        % spikeAmps
        function sa = get.spikeAmps(obj)
            if isfield(obj.dRes, 'spikeAmps') || isprop(obj.dRes, 'spikeAmps')
                sa = obj.dRes.spikeAmps;
            else
                sa = [];
            end
        end

        % spikeClusters/viClu
        function set.spikeClusters(obj, sc)
            obj.spikeClusters = sc;
            obj.nClusters = numel(unique(sc(sc > 0))); %#ok<MCSUP>
        end

        function sc = get.viClu(obj)
            sc = obj.spikeClusters;
        end

        % spikeDelta/delta
        function sd = get.spikeDelta(obj)
            if isfield(obj.sRes, 'spikeDelta')
                sd = obj.sRes.spikeDelta;
            else
                sd = [];
            end
        end
        function sd = get.delta(obj)
            sd = obj.spikeDelta;
        end

        % spikeNeigh/nneigh
        function sn = get.spikeNeigh(obj)
            if isfield(obj.sRes, 'spikeNeigh')
                sn = obj.sRes.spikeNeigh;
            else
                sn = [];
            end
        end
        function sd = get.nneigh(obj)
            sd = obj.spikeNeigh;
        end

        % spikeRho/rho
        function sr = get.spikeRho(obj)
            if isfield(obj.sRes, 'spikeRho')
                sr = obj.sRes.spikeRho;
            else
                sr = [];
            end
        end
        function sr = get.rho(obj)
            sr = obj.spikeRho;
        end

        % spikeFeatures
        function set.spikeFeatures(obj, sf)
            obj.dRes.spikeFeatures = sf;
        end
        function sf = get.spikeFeatures(obj)
            if isfield(obj.dRes, 'spikeFeatures') || isprop(obj.dRes, 'spikeFeatures')
                sf = obj.dRes.spikeFeatures;
            else
                sf = [];
            end
        end

        % spikePositions
        function sf = get.spikePositions(obj)
            if isfield(obj.dRes, 'spikePositions') || isprop(obj.dRes, 'spikePositions')
                sf = obj.dRes.spikePositions;
            else
                sf = [];
            end
        end

        % spikesBySite
        function ss = get.spikesBySite(obj)
            if isfield(obj.dRes, 'spikesBySite') || isprop(obj.dRes, 'spikesBySite')
                ss = obj.dRes.spikesBySite;
            else
                ss = [];
            end
        end

        % spikesBySite2
        function ss = get.spikesBySite2(obj)
            if isfield(obj.dRes, 'spikesBySite2') || isprop(obj.dRes, 'spikesBySite2')
                ss = obj.dRes.spikesBySite2;
            else
                ss = [];
            end
        end

        % spikesBySite3
        function ss = get.spikesBySite3(obj)
            if isfield(obj.dRes, 'spikesBySite3') || isprop(obj.dRes, 'spikesBySite3')
                ss = obj.dRes.spikesBySite3;
            else
                ss = [];
            end
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            ss = obj.dRes.spikeSites;
        end

        % spikeSites2
        function ss = get.spikeSites2(obj)
            if isfield(obj.dRes, 'spikeSites2') || isprop(obj.dRes, 'spikeSites2')
                ss = obj.dRes.spikeSites2;
            else
                ss = [];
            end
        end

        % spikesFilt
        function set.spikesFilt(obj, sf)
            obj.dRes.spikesFilt = sf;
        end
        function sf = get.spikesFilt(obj)
            if isfield(obj.dRes, 'spikesFilt') || isprop(obj.dRes, 'spikesFilt')
                sf = obj.dRes.spikesFilt;
            else
                sf = [];
            end
        end

        % spikesRaw
        function set.spikesRaw(obj, sr)
            obj.dRes.spikesRaw = sr;
        end
        function sr = get.spikesRaw(obj)
            if isfield(obj.dRes, 'spikesRaw') || isprop(obj.dRes, 'spikesRaw')
                sr = obj.dRes.spikesRaw;
            else
                sr = [];
            end
        end

        % spikeTimes
        function ss = get.spikeTimes(obj)
            ss = obj.dRes.spikeTimes;
        end

        % sRes
        function set.sRes(obj, sr)
            if isstruct(sr)
                assert(isfield(sr, 'spikeClusters'), 'spikeClusters must not be missing');
            else
                error('type not recognized: %s', class(sr));
            end

            obj.sRes = sr;
        end

        % tmrWav_clu
        function tm = get.tmrWav_clu(obj)
            tm = obj.meanWfGlobal;
        end

        % unitISIRato/VrIsiRatio_clu
        function ir = get.vrIsiRatio_clu(obj)
            ir = obj.unitISIRatio;
        end

        % unitIsoDist/vrIsoDist_clu
        function id = get.vrIsoDist_clu(obj)
            id = obj.unitIsoDist;
        end

        % unitLRatio/vrLRatio_clu
        function lr = get.vrLRatio_clu(obj)
            lr = obj.unitLRatio;
        end

        % unitPeaks/vrVmin_clu
        function pv = get.vrVmin_clu(obj)
            pv = obj.unitPeaks;
        end

        % unitPeaksRaw/vrVmin_uv_clu
        function pv = get.vrVmin_uv_clu(obj)
            pv = obj.unitPeaksRaw;
        end

        % unitSNR/vrSnr_clu
        function sn = get.vrSnr_clu(obj)
            sn = obj.unitSNR;
        end

        % unitVpp/vrVpp_clu
        function pv = get.vrVpp_clu(obj)
            pv = obj.unitVpp;
        end

        % unitVppRaw/vrVpp_uv_clu
        function pv = get.vrVpp_uv_clu(obj)
            pv = obj.unitVppRaw;
        end

        % vrPos{X,Y}_clu
        function xp = get.vrPosX_clu(obj)
            if ~isempty(obj.clusterCentroids)
                xp = obj.clusterCentroids(:, 1);
            else
                xp = [];
            end
        end
        function yp = get.vrPosY_clu(obj)
            if ~isempty(obj.clusterCentroids)
                yp = obj.clusterCentroids(:, 2);
            else
                yp = [];
            end
        end
    end
end
