classdef Clustering < handle
    %CLUSTERING Model representing clustering of spike data

    %% DETECTION/CLUSTERING RESULTS
    properties (Access=private, Hidden, SetObservable)
        sRes;               % sorting results
        dRes;               % detection results
    end

    properties (SetAccess=private, Hidden, SetObservable)
        hCfg;               % Config object
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
        vnSpk_clu;          % => clusterCounts
        vrIsiRatio_clu;     % => unitISIRatio
        vrIsoDist_clu;      % => unitIsoDist
        vrLRatio_clu;       % => unitLRatio
        vrPosX_clu;         % => clusterCentroids(:, 1)
        vrPosY_clu;         % => clusterCentroids(:, 2)
        vrSnr_clu;          % => unitSNR
        vrVmin_clu;         % => unitPeaks
        vrVmin_uv_clu;      % => unitPeaksRaw
        vrVpp_clu;          % => viSite_clu
        vrVpp_uv_clu;       % => unitVppRaw
        vrVrms_site;        % => siteRMS
    end

    %% NEW-STYLE CLUSTERING PROPERTIES
    properties (SetAccess=private, SetObservable)
        clusterCenters;     % cluster centers
        clusterCentroids;   % centroids of clusters on the probe
        clusterCounts;      % number of spikes per cluster
        clusterNotes;       % notes on clusters
        clusterSites;       % site on which spikes in this cluster most often occur
        history;            % cell array, log of merge/split/delete operations
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
    end

    properties (Hidden)
        meanWfLocal;        % mean filtered waveforms for each cluster
        meanWfGlobal;       % mean filtered waveforms for each cluster over all sites
        meanWfLocalRaw;     % mean raw waveforms for each cluster
        meanWfGlobalRaw;    % mean raw waveforms for each cluster over all sites
        meanWfRawLow;       % mean raw waveforms for each cluster over all sites at a low point on the probe (for drift correction)
        meanWfRawHigh;      % mean raw waveforms for each cluster over all sites at a high point on the probe (for drift correction)
        simScore;           % waveform-based cluster similarity scores
        tmrWav_clu;
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

    % computed from other values
    properties (Dependent, Transient)
        nClusters;          % number of clusters
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
        function obj = Clustering(sRes, dRes, hCfg)
            obj.sRes = sRes;
            obj.dRes = dRes;
            obj.hCfg = hCfg;

            % these fields are mutable so we need to store copies in obj
            obj.spikeClusters = obj.initialClustering;
            if isfield(sRes, 'clusterCenters')
                obj.clusterCenters = sRes.clusterCenters;
            else
                obj.clusterCenters = [];
            end
            obj.clusterCentroids = [];
            obj.history = cell(0, 2);

            obj.clearNotes();
            obj.refresh(true);
            obj.commit('initial commit');
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
                fprintf('Computing self correlation\n\t');
                t1 = tic;

                selfSim = zeros(1, obj.nClusters);
                for iCluster = 1:obj.nClusters
                    selfSim(iCluster) = scKern(obj, iCluster);
                    fprintf('.');
                end

                fprintf('\n\ttook %0.1fs\n', toc(t1));
            else
                selfSim = scKern(obj, iCluster);
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
            %AUTOMERGEBYSIM Automatically merge clusters by similarity score
            simScore_ = obj.simScore;
            nClusters_ = size(simScore_, 2);

            % Identify clusters to remove, update and same (no change), disjoint sets
            simScore_(tril(true(nClusters_))) = 0; % ignore bottom half
            [scoresMax, mapTo] = max(simScore_);
            keepMe_ = scoresMax < obj.hCfg.maxWavCor; % keep clusters whose max similarity to another cluster is less than our threshold

            if all(keepMe_)
                nMerged = 0;
                return;
            end

            minScore = min(scoresMax(~keepMe_));
            keepMe = find(keepMe_);
            mapTo(keepMe_) = keepMe;
            keepMe = setdiff(keepMe, mapTo(~keepMe_));
            removeMe = setdiff(1:nClusters_, mapTo);
            updateMe = setdiff(setdiff(1:nClusters_, keepMe), removeMe);

            % update cluster number
            % try
            %     obj.clusterCenters(removeMe) = [];
            % catch ME
            %     rethrow(ME);
            % end

            %hClust = S_clu_map_index_(hClust, mapTo); %index mapped
            good = obj.spikeClusters > 0;
            mapTo = int32(mapTo);
            obj.spikeClusters(good) = mapTo(obj.spikeClusters(good)); %translate cluster number
            obj.refresh(false); % empty clusters removed later

            obj.rmRefracSpikes(); % remove refrac spikes

            % update cluster waveforms and distance
            obj.computeMeanWaveforms(updateMe);
            obj.computeWaveformSim();
            obj.removeEmptyClusters();

            nMerged = nClusters_ - obj.nClusters;
            if obj.hCfg.verbose
                fprintf('\tnClusters: %d->%d (%d merged, min score: %0.4f)\n', nClusters_, obj.nClusters, nMerged, minScore);
            end
        end

        function removeEmptyClusters(obj)
            %REMOVEEMPTYCLUSTERS Find and remove empty clusters
            keepClusters = obj.clusterCounts > 0;
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

        function subsetFields(obj, keepMe) % WIP
            %SELECT Subset all data fields, taking only those indices we want to keep
            % subset vector fields
            if ~isempty(obj.clusterCenters) % ?
                obj.clusterCenters = obj.clusterCenters(keepMe);
            end
            if ~isempty(obj.clusterCounts) % ?
                obj.clusterCounts = obj.clusterCounts(keepMe);
            end
            if ~isempty(obj.clusterSites) % ?
                obj.clusterSites = obj.clusterSites(keepMe);
            end
            if ~isempty(obj.nSitesOverThresh)
                obj.nSitesOverThresh = obj.nSitesOverThresh(keepMe);
            end
            if ~isempty(obj.unitISIRatio)
                obj.unitISIRatio = obj.unitISIRatio(keepMe);
            end
            if ~isempty(obj.unitIsoDist)
                obj.unitIsoDist = obj.unitIsoDist(keepMe);
            end
            if ~isempty(obj.unitLRatio)
                obj.unitLRatio = obj.unitLRatio(keepMe);
            end
            if ~isempty(obj.unitPeaks)
                obj.unitPeaks = obj.unitPeaks(keepMe);
            end
            if ~isempty(obj.unitPeaksRaw)
                obj.unitPeaksRaw = obj.unitPeaksRaw(keepMe);
            end
            if ~isempty(obj.unitPeakSites)
                obj.unitPeakSites = obj.unitPeakSites(keepMe);
            end
            if ~isempty(obj.unitSNR)
                obj.unitSNR = obj.unitSNR(keepMe);
            end
            if ~isempty(obj.unitVpp)
                obj.unitVpp = obj.unitVpp(keepMe);
            end
            if ~isempty(obj.unitVppRaw)
                obj.unitVppRaw = obj.unitVppRaw(keepMe);
            end

            % subset matrix fields
            if ~isempty(obj.clusterCentroids) % nClusters x 2 %?
                obj.clusterCentroids = obj.clusterCentroids(keepMe, :);
            end
            if ~isempty(obj.simScore) % nClusters x nClusters
                obj.simScore = obj.simScore(keepMe, keepMe);
            end

            % subset tensor fields
            if ~isempty(obj.meanWfGlobal)
                obj.meanWfGlobal = obj.meanWfGlobal(:, :, keepMe);
            end
            if ~isempty(obj.meanWfGlobalRaw)
                obj.meanWfGlobalRaw = obj.meanWfGlobalRaw(:, :, keepMe);
            end
            if ~isempty(obj.meanWfLocal)
                obj.meanWfLocal = obj.meanWfLocal(:, :, keepMe);
            end
            if ~isempty(obj.meanWfLocalRaw)
                obj.meanWfLocalRaw = obj.meanWfLocalRaw(:, :, keepMe);
            end
            if ~isempty(obj.meanWfRawHigh)
                obj.meanWfRawHigh = obj.meanWfRawHigh(:, :, keepMe);
            end
            if ~isempty(obj.meanWfRawLow)
                obj.meanWfRawLow = obj.meanWfRawLow(:, :, keepMe);
            end

            % subset cell fields
            if ~isempty(obj.clusterNotes)
                obj.clusterNotes = obj.clusterNotes(keepMe);
            end
            if ~isempty(obj.spikesByCluster)
                obj.spikesByCluster = obj.spikesByCluster(keepMe);
            end
        end
    end

    %% USER METHODS
    methods
        function autoMerge(obj, doAssign)
            %AUTOMERGE Automatically merge clusters
            if nargin < 4
                doAssign = true;
            end

        %     if doAssign
        %         S_clu = postCluster_(S_clu, P);
        %     end

            obj.refresh(true);
            obj.orderClusters('clusterSites');
            obj.clearNotes();

            obj.rmOutlierSpikes();
            obj.doWaveformMerge();
            obj.refresh(true);
            obj.orderClusters('clusterSites');
            obj.updateWaveforms();

            obj.computeCentroids();
            obj.clearNotes();
            obj.computeQualityScores();
            obj.commit('autoMerge');
        end

        function clearNotes(obj)
            %CLEARNOTES Remove all cluster notes
            obj.clusterNotes = cell(obj.nClusters, 1);
        end

        function commit(obj, msg)
            %COMMIT Commit a modification of clustering to history log
            if nargin < 2
                msg = '';
            end

            if isempty(obj.history)
                iDiffs = find(obj.spikeClusters ~= obj.initialClustering);
                sDiffs = obj.spikeClusters(iDiffs);
                obj.history(1, :) = {msg, [iDiffs'; sDiffs']};
            else
                % check for consistency before committing
                ic = obj.selfConsistent();
                if ~isempty(ic)
                    wmsg = strjoin(ic, '\n\t');
                    warning('Cluster data inconsistent after previous operations:\n\t%s', wmsg);
                    obj.editSeek(obj.nEdits);
                    return;
                end

                % replay from the beginning
                spikeClusters_ = obj.initialClustering;

                for j = 1:size(obj.history, 1)
                    diffs = obj.history{j, 2};
                    if isempty(diffs)
                        continue;
                    end

                    iDiffs = diffs(1, :);
                    sDiffs = diffs(2, :);
                    spikeClusters_(iDiffs) = sDiffs;
                end % now caught up to the current state, can get diff
                
                iDiffs = find(obj.spikeClusters ~= spikeClusters_);
                sDiffs = obj.spikeClusters(iDiffs);
                obj.history(end+1, :) = {msg, [iDiffs'; sDiffs']};
            end
        end

        function computeCentroids(obj, updateMe)
            %COMPUTEPOSITIONS determine cluster position from spike position
            if nargin < 2 || isempty(obj.centroids)
                updateMe = [];
            end

            if isempty(updateMe)
                centroids_ = zeros(obj.nClusters, 2);
                clusters_ = 1:obj.nClusters;
            else % selective update
                centroids_ = obj.clusterCentroids;
                clusters_ = updateMe(:)';
            end

            featureSites = 1:(1+obj.hCfg.nSiteDir*2 - obj.hCfg.nSitesExcl);

            for iCluster = clusters_
                [clusterSpikes, neighbors] = subsampleCenteredSpikes(obj, iCluster);
                if isempty(clusterSpikes)
                    continue;
                end

                neighbors = neighbors(1:end-obj.hCfg.nSitesExcl);
                featureWeights = squeeze(obj.spikeFeatures(featureSites, 1, clusterSpikes));
                neighborLoc = single(obj.hCfg.siteLoc(neighbors, :)); % position on probe

                centroids_(iCluster, 1) = median(getWeightedLoc(featureWeights, neighborLoc(:, 1)));
                centroids_(iCluster, 2) = median(getWeightedLoc(featureWeights, neighborLoc(:, 2)));
            end

            obj.clusterCentroids = centroids_;
        end

        function computeMeanWaveforms(obj, updateMe, useRaw)
            %COMPUTEMEANWAVEFORMS Compute the mean waveform for each cluster
            if nargin < 2
                updateMe = [];
            end
            if nargin < 3
                useRaw = true;
            end

            results = compMeanWf(obj, updateMe, useRaw);
            
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

        function computeWaveformSim(obj, updateMe_)
            %WAVEFORMSIM Compute waveform-based similarity scores for all clusters
            if nargin < 3 || isempty(obj.simScore)
                updateMe_ = [];
            end

            if isempty(obj.meanWfGlobal)
                return;
            end

            obj.hCfg.useGPU = false; % disable GPU for this
            obj.simScore = compWfSim(obj, updateMe_);
            obj.hCfg.useGPU = true; % disable GPU for this
        end

        function computeQualityScores(obj, updateMe)
            %COMPUTEQUALITYSCORES Get cluster quality scores
            if nargin < 2
                updateMe = [];
            end
            scores = qualScores(obj, updateMe);
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

        function doWaveformMerge(obj, maxWavCor)
            %DOWAVEFORMMERGE Merge clusters by waveform-based similarity scores
            if nargin == 2 && ~isempty(maxWavCor)
                mwcOld = obj.hCfg.maxWavCor;
                obj.hCfg.maxWavCor = mwcOld;
            end

            obj.computeMeanWaveforms();
            obj.computeWaveformSim();

            for iRepeat = 1:obj.hCfg.nPassesMerge % single-pass vs dual-pass correction
                nMerged = obj.mergeBySim();
                if nMerged < 1
                    break;
                end
            end

            if nargin == 2 && ~isempty(maxWavCor) % restore old maxWavCor
                obj.hCfg.maxWavCor = mwcOld;
            end
        end

        function editSeek(obj, seekTo)
            %EDITSEEK Seek back and forth to position in history
            %   TODO: this can cause branching issues; make non-tip
            %   positions read-only
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
                diffs = obj.history{j, 2};
                if isempty(diffs)
                    continue;
                end

                iDiffs = diffs(1, :);
                sDiffs = diffs(2, :);
                spikeClusters_(iDiffs) = sDiffs;
            end

            obj.spikeClusters = spikeClusters_;

            obj.refresh(true);
            if ~isempty(obj.meanWfGlobal) % only compute mean waveforms if we already had them
                obj.updateWaveforms();
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

        function inconsistencies = selfConsistent(obj)
            %SELFCONSISTENT Check all fields have the correct sizes
            inconsistencies = {};
            % vector fields
            if ~isempty(obj.clusterCenters) && numel(obj.clusterCenters) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterCenters: expected %d, actual %d', obj.nClusters, numel(obj.clusterCenters));
            end
            if ~isempty(obj.clusterCounts) && numel(obj.clusterCounts) ~= obj.nClusters
                inconsistencies{end+1} = sprintf('clusterCounts: expected %d, actual %d', obj.nClusters, numel(obj.clusterCounts));
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
                inconsistencies{end+1} = sprintf('simScore: expected %dx%d, actual %d', obj.nClusters, obj.nClusters, size(obj.simScore, 1), size(obj.simScore, 2));
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
        end

        function orderClusters(obj, by)
            %ORDERCLUSTERS Arrange cluster ID numbers by some criterion

            if nargin < 2 || isempty(by)
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
            obj.spikesByCluster = obj.spikesByCluster(argsort); % redundant with refresh?
            obj.clusterSites = obj.clusterSites(argsort);       % redundant with refresh?
            obj.clusterCounts = obj.clusterCounts(argsort);     % redundant with refresh?
            obj.clusterNotes = obj.clusterNotes(argsort);
            if ~isempty(obj.clusterCentroids)
                obj.clusterCentroids = obj.clusterCentroids(argsort, :);
            end

            obj.refresh(true);
        end

        function refresh(obj, doRemoveEmpty)
            %REFRESH Recount and store spikes by cluster, optionally removing empty clusters
            obj.spikesByCluster = arrayfun(@(iC) find(obj.spikeClusters == iC), (1:obj.nClusters)', 'UniformOutput', false);
            obj.clusterCounts = cellfun(@numel, obj.spikesByCluster);
            if ~isempty(obj.spikeSites)
                obj.clusterSites = double(arrayfun(@(iC) mode(obj.spikeSites(obj.spikesByCluster{iC})), 1:obj.nClusters));
            else
                obj.clusterSites = [];
            end

            if doRemoveEmpty
                obj.removeEmptyClusters();
            end
        end

        function rmOutlierSpikes(obj)
            %RMOUTLIERSPIKES Mahalanobis-distance based outlier removal
            if obj.hCfg.outlierThresh == 0
                return;
            end

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
                    iDist = jrclust.utils.madScore(log(mahal(iFeatures, iFeatures)));
                catch
                    continue;
                end

                exclude = (iDist > obj.hCfg.outlierThresh);
                if any(exclude)
                    obj.spikesByCluster{iCluster} = iSpikes(~exclude);
                    obj.spikeClusters(iSpikes(exclude)) = 0; % classify as noise
                    obj.clusterCounts(iCluster) = numel(obj.spikesByCluster{iCluster});
                end

                fprintf('.');
            end
        end

        function nRemoved = rmRefracSpikes(obj, iCluster)
            % remove refractory spikes
            if nargin == 1 % recurse for each cluster
                nRemoved = 0;
                for iCluster_ = 1:obj.nClusters
                    nRemoved_ = obj.rmRefracSpikes(iCluster_);
                    nRemoved = nRemoved + nRemoved_;
                end

                return;
            else
                nRemoved = 0;
                nSkip_refrac = 4;

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

                    keepMe(iKeepMe(inRefrac(1:nSkip_refrac:end))) = false;
                end

                nRemoved = sum(~keepMe);
                nTotal = numel(keepMe);
                obj.spikeClusters(clusterSpikes_(~keepMe)) = 0; % assign to noise cluster

                obj.spikesByCluster{iCluster} = clusterSpikes_(keepMe);
                obj.clusterCounts(iCluster) = sum(keepMe);
            end

            if obj.hCfg.verbose
                fprintf('Cluster %d: removed %d/%d (%0.2f%%) duplicate spikes\n', iCluster, nRemoved, nTotal, 100*nRemoved/nTotal);
            end
        end

        function updateWaveforms(obj)
            obj.computeMeanWaveforms();
            obj.computeWaveformSim();
            obj.simScore = jrclust.utils.setDiag(obj.mrWavCor, obj.computeSelfSim());
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

        % clusterCounts/vnSpk_clu
        function cc = get.vnSpk_clu(obj)
            cc = obj.clusterCounts;
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
        function nc = get.nClusters(obj)
            nc = double(max(obj.spikeClusters));
        end
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
        function sf = get.spikesFilt(obj)
            if isfield(obj.dRes, 'spikesFilt') || isprop(obj.dRes, 'spikesFilt')
                sf = obj.dRes.spikesFilt;
            else
                sf = [];
            end
        end

        % spikesRaw
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
