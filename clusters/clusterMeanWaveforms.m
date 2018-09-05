%--------------------------------------------------------------------------
% 10/10/17 JJJ: moved spikeWaveforms and spikeTraces internally
function S_clu = clusterMeanWaveforms(S_clu, clustersToUpdate, fSkipRaw)
    % average cluster waveforms and determine the center
    % only use the centered spikes
    % viClu_copy: waveforms not changing

    if nargin < 2
        clustersToUpdate = [];
    end
    if nargin < 3
        fSkipRaw = 0;
    end

    S0 = get(0, 'UserData');
    P = S0.P;

    P.fMeanSubt = 0;
    fVerbose = isempty(clustersToUpdate);

    if fVerbose
        fprintf('Calculating cluster mean waveform.\n\t');
        t1 = tic;
    end

    if isfield(S_clu, 'nClusters')
        nClusters = S_clu.nClusters;
    else
        nClusters = max(S_clu.spikeClusters);
    end

    nSamples = S0.waveformDims(1);
    nSitesProbe = numel(P.chanMap); % number of sites on the entire probe
    nSitesSpike = S0.waveformDims(2); % number of sites per event

    % Prepare cluster loop
    trWav_spk_clu = zeros(nSamples, nSitesSpike, nClusters, 'single');
    tmrWav_spk_clu = zeros(nSamples, nSitesProbe, nClusters, 'single');

    if ~fSkipRaw
        rawTraces = []; % clear memory
        rawTraces = getSpikeWaveforms(P, 1);

        nSamplesRaw = S0.traceDims(1);
        trWav_raw_clu = zeros(nSamplesRaw, nSitesSpike, nClusters, 'single');

        tmrWav_raw_clu = zeros(nSamplesRaw, nSitesProbe, nClusters, 'single');
        tmrWav_raw_lo_clu = zeros(nSamplesRaw, nSitesProbe, nClusters, 'single');
        tmrWav_raw_hi_clu = zeros(nSamplesRaw, nSitesProbe, nClusters, 'single');
    else
        trWav_raw_clu = [];
        tmrWav_raw_clu = [];
        tmrWav_raw_lo_clu = [];
        tmrWav_raw_hi_clu = [];
    end

    if ~isempty(clustersToUpdate) % specific clusters we want to update
        updateCluster = ismember(1:nClusters, clustersToUpdate)'; % clusters
        nClustersPrevious = size(S_clu.trWav_spk_clu, 3);
        updateCluster((1:nClusters) > nClustersPrevious) = 1;

        % load waveform data from S_clu to alter
        tmrWav_spk_clu = S_clu.tmrWav_spk_clu;
        tmrWav_raw_clu = S_clu.tmrWav_raw_clu;

        trWav_spk_clu = S_clu.trWav_spk_clu;
        trWav_raw_clu = S_clu.trWav_raw_clu;

        tmrWav_raw_lo_clu = S_clu.tmrWav_raw_lo_clu;
        tmrWav_raw_hi_clu = S_clu.tmrWav_raw_hi_clu;
    else % update everything, don't bother loading data from S_clu
        updateCluster = true(nClusters, 1);
    end

    filteredTraces = getSpikeWaveforms(P, 0);

    for iCluster = 1:nClusters
        if updateCluster(iCluster) % compute mean filtered waveforms
            [meanTraces, iClusterSites] = meanWaveforms(S_clu, filteredTraces, iCluster, S0);
            if isempty(meanTraces)
                continue;
            end

            % cluster-mean filtered waveforms just in the neighborhood of the center site
            trWav_spk_clu(:, :, iCluster) = bit2uV_(meanTraces, P)
            % cluster-mean filtered waveforms across all sites (0 outside neighborhood)
            tmrWav_spk_clu(:, iClusterSites, iCluster) = trWav_spk_clu(:, :, iCluster);

            if ~fSkipRaw % also compute mean raw waveforms
                [meanTraces, iClusterSites, mrWav_lo_clu1, mrWav_hi_clu1] = meanWaveforms(S_clu, rawTraces, iCluster, S0);
                if isempty(meanTraces)
                    continue;
                end

                % cluster-mean raw traces just in the neighborhood of the center site
                trWav_raw_clu(:, :, iCluster) = meanTraces*P.uV_per_bit;
                % cluster-mean raw traces across all sites (0 outside neighborhood)
                tmrWav_raw_clu(:, iClusterSites, iCluster) = trWav_raw_clu(:, :, iCluster);

                if isempty(mrWav_lo_clu1) || isempty(mrWav_hi_clu1)
                    tmrWav_raw_lo_clu(:, iClusterSites, iCluster) = zeros(nSamplesRaw, numel(iClusterSites));
                    tmrWav_raw_hi_clu(:, iClusterSites, iCluster) = zeros(nSamplesRaw, numel(iClusterSites));
                else
                    tmrWav_raw_lo_clu(:, iClusterSites, iCluster) = mrWav_lo_clu1*P.uV_per_bit;
                    tmrWav_raw_hi_clu(:, iClusterSites, iCluster) = mrWav_hi_clu1*P.uV_per_bit;
                end
            end
        end
        if fVerbose
            fprintf('.');
        end
    end % cluster

    tmrWav_clu = tmrWav_spk_clu;

    % get min waveform values per cluster and on which site they occur
    [vrVmin_clu, viSite_min_clu] = min(squeeze(min(trWav_spk_clu)), [], 1);
    vrVmin_clu = abs(vrVmin_clu(:)); % min waveform values per cluster
    viSite_min_clu = viSite_min_clu(:); % site on which min occurs

    S_clu = struct_add_(S_clu, vrVmin_clu, viSite_min_clu, ...
        trWav_spk_clu, tmrWav_spk_clu, trWav_raw_clu, tmrWav_raw_clu, tmrWav_clu, ...
        tmrWav_raw_lo_clu, tmrWav_raw_hi_clu);

    if fVerbose
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    end
end % function

%--------------------------------------------------------------------------
% 10/22/17 JJJ
function [meanTraces, clusterSites, mrWav_lo_clu1, mrWav_hi_clu1] = meanWaveforms(S_clu, traces, cluster, S0)
    % compute the mean-subtracted average waveform of CLUSTER

    if nargin < 4
        S0 = get(0, 'UserData');
    end

    nSamples_max = 1000;
    fUseCenterSpk = getOr(S0.P, 'fUseCenterSpk', 0); % set to zero to use all spikes
    fDrift_merge = getOr(S0.P, 'fDrift_merge', 0);

    meanTraces = [];
    mrWav_lo_clu1 = [];
    mrWav_hi_clu1 = [];

    centerSite = S_clu.clusterSites(cluster);
    clusterSites = S0.P.miSites(:, centerSite); % neighborhood of cluster center site
    clusterSpikes = S_clu.spikesByCluster{cluster};
    clusterSpikeSites = S0.spikeSites(clusterSpikes); % center sites of spikes in cluster

    if fUseCenterSpk % only take spikes occurring on the center site
        isCentered = (clusterSpikeSites == centerSite);
        clusterSpikes = clusterSpikes(isCentered);
        clusterSpikeSites = clusterSpikeSites(isCentered);
    end

    if isempty(clusterSpikes)
        return;
    end

    if ~fDrift_merge
        spikes = spikesNearMidpoint(clusterSpikes, S0.spikeTimes, 1/S0.P.nTime_clu);
        meanTraces = mean(single(traces(:, :, spikes)), 3);
        meanTraces = meanSubtract(meanTraces); %122717 JJJ

        return;
    end

    % TODO: pick up here -- acl
    if ~isfield(S0, 'mrPos_spk')
        S0.mrPos_spk = spk_pos_(S0);
    end
    vrPosY_spk1 = S0.mrPos_spk(clusterSpikes, 2); % position based quantile
    vrYLim = quantile(vrPosY_spk1, [0,1,2,3]/3);
    [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(2:3), nSamples_max, clusterSpikeSites);
    meanTraces = nanmean_int16_(traces(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P); % * S0.P.uV_per_bit;

    if nargout > 2
        [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(1:2), nSamples_max, clusterSpikeSites);
        mrWav_lo_clu1 = nanmean_int16_(traces(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P);

        [viSpk_clu_, clusterSites_] = spk_select_pos_(clusterSpikes, vrPosY_spk1, vrYLim(3:4), nSamples_max, clusterSpikeSites);
        mrWav_hi_clu1 = nanmean_int16_(traces(:,:,viSpk_clu_), 3, fUseCenterSpk, centerSite, clusterSites_, S0.P);
    end
end % function
