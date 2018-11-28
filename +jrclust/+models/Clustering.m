classdef Clustering < handle
    %CLUSTERING Model representing clustering of spike data

    %% DETECTION/CLUSTERING RESULTS
    properties (Access=private, Hidden, SetObservable)
        sRes;               % sorting results
        dRes;               % detection results
    end

    %% OLD-STYLE PROPERTIES, publicly gettable (will be deprecated after a grace period)
    properties (SetObservable, Dependent, Hidden, Transient)
        cviSpk_clu;         % => spikesByCluster
        csNote_clu;         % => clusterNotes
        nClu;               % => nClusters
        tmrWav_spk_clu;     % => meanWfGlobal
        tmrWav_raw_clu;     % => meanWfGlobalRaw
        trWav_spk_clu;      % => meanWfLocal
        trWav_raw_clu;      % => meanWfLocalRaw
        tmrWav_raw_hi_clu;  % => meanWfRawHigh
        tmrWav_raw_lo_clu;  % => meanWfRawLow
        viClu;              % => spikeClusters
        viSite_clu;         % => clusterSites
        viSite_min_clu;     % => minSites
        vnSpk_clu;          % => clusterCounts
        vrPosX_clu;         % => centroids(:, 1)
        vrPosY_clu;         % => centroids(:, 2)
        vrVmin_clu;         % => peakVoltages
    end

    %% NEW-STYLE CLUSTERING PROPERTIES
    properties (Dependent, Transient)
        initialClustering;  % initial assignment of spikes to cluster
        nClusters;          % number of clusters
    end

    properties (SetAccess=private)
        clusterCounts;      % number of spikes per cluster
        clusterSites;       % site on which spikes in this cluster most often occur
        clusterNotes;       % notes on clusters
        spikeClusters;      % individual spike assignments
        spikesByCluster;    % cell array of spike indices per cluster
    end

    properties (Hidden)
        peakVoltages; 
        minSites;

        meanWfLocal;
        meanWfGlobal;
        meanWfLocalRaw;
        meanWfGlobalRaw;
        meanWfRawLow;
        meanWfRawHigh;

        mrWavCor;
        tmrWav_clu;
    end

    properties (SetObservable)
        centroids;          % centroids of clusters on the probe
    end

    %% DETECTION RESULTS
    properties
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
        function obj = Clustering(sRes, dRes)
            obj.sRes = sRes;
            obj.dRes = dRes;

            obj.spikeClusters = obj.initialClustering;
            obj.clearNotes();
            obj.refresh(true);
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = getClusterMean(obj, spikeWindows, iCluster, hCfg)
            %GETCLUSTERMEAN Get mean cluster waveform, optionally low and high

            function [spikesOut, sitesOut] = selectSpikesInBounds(spikesIn, sitesIn, yPos, yLims)
                %SELECTSPIKESINBOUNDS Select spikes by whether their y-positions fall within given limits
                nSamplesMax = 1000;
                inBounds = yPos >= yLims(1) & yPos < yLims(2);

                if ~any(inBounds)
                    spikesOut = spikesIn;
                    sitesOut = sitesIn;
                    return;
                end

                spikesOut = jrclust.utils.subsample(spikesIn(inBounds), nSamplesMax);
                sitesOut = jrclust.utils.subsample(sitesIn(inBounds), nSamplesMax);
            end

            function mrWav_clu1 = nanmeanInt16(spikeWindows, iSite, sites, hCfg)
                iSiteNeighbors = hCfg.siteNeighbors(:, iSite);
                trWav = nan([size(spikeWindows, 1), numel(iSiteNeighbors), numel(sites)], 'single');
                uniqueSites = unique(sites);
                nUniqueSites = numel(uniqueSites);
                uniqueNeighbors = hCfg.siteNeighbors(:, uniqueSites);

                for jSite = 1:nUniqueSites
                    iSiteUnique = uniqueSites(jSite);
                    viSpk_ = find(sites == iSiteUnique);

                    [~, viSite1a_, viSite1b_] = intersect(iSiteNeighbors, uniqueNeighbors(:, jSite));
                    if isempty(viSite1a_)
                        continue;
                    end

                    trWav(:, viSite1a_, viSpk_) = spikeWindows(:, viSite1b_, viSpk_);
                end

                mrWav_clu1 = nanmean(trWav, 3);
                mrWav_clu1 = jrclust.utils.meanSubtract(mrWav_clu1); %122717 JJJ
            end

            [clusterMean, clusterMeanLow, clusterMeanHigh] = deal([]);
            iSite = obj.clusterSites(iCluster);
            siteNeighbors = hCfg.siteNeighbors(:, iSite);

            clusterSpikes = obj.spikesByCluster{iCluster};
            clusterSites_ = obj.spikeSites(clusterSpikes);

            if isempty(clusterSpikes)
                return;
            end

            if ~hCfg.fDrift_merge || isempty(obj.spikePositions)
                middlemost = spk_select_mid_(clusterSpikes, obj.spikeTimes, hCfg.nTime_clu);
                clusterMean = mean(single(spikeWindows(:, :, middlemost)), 3);
                clusterMean = jrclust.utils.meanSubtract(clusterMean);
                return;
            end

            yPos = obj.spikePositions(clusterSpikes, 2); % position based quantile
            yLims = quantile(yPos, [0, 1, 2, 3]/3);

            [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(2:3));
            clusterMean = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hCfg); % * hCfg.uV_per_bit;

            if nargout > 2
                [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(1:2));
                clusterMeanLow = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hCfg);

                [selectedSpikes, selectedSites] = selectSpikesInBounds(clusterSpikes, clusterSites_, yPos, yLims(3:4));
                clusterMeanHigh = nanmeanInt16(spikeWindows(:, :, selectedSpikes), iSite, selectedSites, hCfg);
            end
        end

        function removeEmptyClusters(obj)
            %REMOVEEMPTYCLUSTERS Find and remove empty clusters
            keepClusters = obj.clusterCounts > 0;
            if all(keepClusters)
                return;
            end

            % subset all fields indexed by cluster
            % obj.subsetFields(obj, keepClusters);

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

        function subsetFields(obj, keepClusters) % WIP
            %SELECT Subset all data fields, taking only those we're interested in
            % excl vnSpk_clu, viSite_clu, vrPosX_clu, vrPosY_clu

            csNames = fieldnames(obj);
            if isempty(csNames)
                return;
            end

            % subset vector fields
            viMatch_v = cellfun(@(vi) ~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_v), keepClusters);

            % subset tensor fields
            viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_t), keepClusters, 3);

            % subset cell fields
            viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
            obj = struct_select_(obj, csNames(viMatch_c), keepClusters);

            % remap mrWavCor
            if isprop(obj, 'mrWavCor')
        %         obj.mrWavCor = S_clu_wavcor_remap_(obj, viKeep_clu);
                obj.mrWavCor = obj.mrWavCor(keepClusters, keepClusters);
            end

            % remap mrSim_clu
            if isprop(obj, 'mrSim_clu')
                obj.mrSim_clu = obj.mrSim_clu(keepClusters, keepClusters);
            end
        end
    end

    %% USER METHODS
    methods
        function clearNotes(obj)
            %CLEARNOTES Remove all cluster notes
            obj.clusterNotes = cell(obj.nClusters, 1);
        end

        function computeMeanWaveforms(obj, hCfg, updateMe, useRaw)
            %COMPUTEMEANWAVEFORMS Compute the mean waveform for each cluster
            if isempty(obj.spikesFilt)
                return;
            end

            if nargin < 3
                updateMe = [];
            end
            if nargin < 4
                useRaw = true;
            end
            useRaw = useRaw && ~isempty(obj.spikesRaw);
            
            verbose = hCfg.verbose && isempty(updateMe);
            if verbose
                fprintf('Calculating cluster mean waveform.\n\t');
                t = tic;
            end

            nSites = numel(hCfg.siteMap);
            [nSamples, nSitesEvt, ~] = size(obj.spikesFilt);

            meanWfLocal_ = zeros(nSamples, nSitesEvt, obj.nClusters, 'single');
            meanWfGlobal_ = zeros(nSamples, nSites, obj.nClusters, 'single');

            if useRaw
                nSamplesRaw = size(obj.spikesRaw, 1);
                meanWfLocalRaw_ = zeros(nSamplesRaw, nSitesEvt, obj.nClusters, 'single');
                meanWfGlobalRaw_ = zeros(nSamplesRaw, nSites, obj.nClusters, 'single');
                meanWfRawLow_ = zeros(nSamplesRaw, nSites, obj.nClusters, 'single');
                meanWfRawHigh_ = zeros(nSamplesRaw, nSites, obj.nClusters, 'single');
            else
                [meanWfLocalRaw_, meanWfGlobalRaw_, meanWfRawLow_, meanWfRawHigh_] = deal([]);
            end

            % we have specific clusters to update
            if ~isempty(updateMe) && ~isempty(obj.meanWfLocal)
                % visit all clusters explicitly requested or not previously seen
                visitMe = false(obj.nClusters, 1);
                visitMe(updateMe) = true;
                visitMe((1:obj.nClusters) > size(obj.meanWfLocal, 3)) = true;

                meanWfLocal_ = obj.meanWfLocal;
                meanWfGlobal_ = obj.meanWfGlobal;

                meanWfLocalRaw_ = obj.meanWfLocalRaw;
                meanWfGlobalRaw_ = obj.meanWfGlobalRaw;

                meanWfRawLow_ = obj.tmrWav_raw_lo_clu;
                meanWfRawHigh_ = obj.tmrWav_raw_hi_clu;
            else % nothing specific requested, update all
                visitMe = true(obj.nClusters, 1);
            end

            for iCluster = 1:obj.nClusters
                if visitMe(iCluster)
                    [clusterMean, siteNeighbors] = getClusterMean(obj, obj.spikesFilt, iCluster, hCfg);

                    if isempty(clusterMean)
                        continue;
                    end

                    clusterMean = jrclust.utils.bit2uV(clusterMean, hCfg);
                    meanWfLocal_(:, :, iCluster) = clusterMean;
                    meanWfGlobal_(:, siteNeighbors, iCluster) = clusterMean;
                end

                if verbose
                    fprintf('.');
                end
            end

            % Compute spkraw
            if useRaw
                for iCluster = 1:obj.nClusters
                    if visitMe(iCluster)
                        [clusterMean, siteNeighbors, clusterMeanLow, clusterMeanHigh] = getClusterMean(obj, obj.spikesRaw, iCluster, hCfg);
                        if isempty(clusterMean)
                            continue;
                        end

                        clusterMean = jrclust.utils.meanSubtract(clusterMean)*hCfg.bitScaling;
                        meanWfGlobalRaw_(:, siteNeighbors, iCluster) = clusterMean;
                        meanWfLocalRaw_(:, :, iCluster) = clusterMean;

                        if isempty(clusterMeanLow) || isempty(clusterMeanHigh)
                            meanWfRawLow_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                            meanWfRawHigh_(:, siteNeighbors, iCluster) = zeros(nSamplesRaw, numel(siteNeighbors));
                        else
                            meanWfRawLow_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanLow)*hCfg.bitScaling;
                            meanWfRawHigh_(:,siteNeighbors,iCluster) = jrclust.utils.meanSubtract(clusterMeanHigh)*hCfg.bitScaling;
                        end
                    end

                    if verbose
                        fprintf('.');
                    end
                end
            end

            obj.tmrWav_clu = meanWfGlobal_; % meanSubtract after or before?

            % measure waveforms
            [peakVoltages_, minSites_] = min(permute(min(meanWfLocal_), [2, 3, 1]), [], 1);
            obj.peakVoltages = abs(peakVoltages_(:));
            obj.minSites = minSites_(:);

            % collect computed values
            obj.meanWfLocal = meanWfLocal_;
            obj.meanWfGlobal = meanWfGlobal_;
            obj.meanWfLocalRaw = meanWfLocalRaw_;
            obj.meanWfGlobalRaw = meanWfGlobalRaw_;
            obj.meanWfRawLow = meanWfRawLow_;
            obj.meanWfRawHigh = meanWfRawHigh_;

            if verbose
                fprintf('\n\ttook %0.1fs\n', toc(t));
            end
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

        function orderClusters(obj, by)
            %ORDERCLUSTERS Arrange cluster ID numbers by some criterion

            if nargin < 2 || isempty(by)
                by = 'clusterSites';
            end

            if strcmpi(by, 'Y + X') && ~isempty(obj.centroids)
                [~, argsort] = sort(sum(obj.centroids, 2), 'ascend');
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
            if ~isempty(obj.centroids)
                obj.centroids = obj.centroids(argsort, :);
            end

            obj.refresh();
        end
    end

    %% GETTERS/SETTERS
    methods
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

        % initialClustering
        function ic = get.initialClustering(obj)
            if ~isempty(obj.sRes)
                ic = obj.sRes.spikeClusters;
            else
                ic = [];
            end
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

        % minSites/viSite_min_clu
        function ms = get.viSite_min_clu(obj)
            ms = obj.minSites;
        end

        % nClusters/nClu
        function nc = get.nClusters(obj)
            nc = double(max(obj.spikeClusters));
        end
        function nc = get.nClu(obj)
            nc = obj.nClusters;
        end

        % peakVoltages/vrVmin_clu
        function pv = get.vrVmin_clu(obj)
            pv = obj.peakVoltages;
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
    end
end
