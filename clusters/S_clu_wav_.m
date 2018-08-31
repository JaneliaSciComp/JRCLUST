%--------------------------------------------------------------------------
% 10/10/17 JJJ: moved spikeWaveforms and spikeTraces internally
function S_clu = S_clu_wav_(S_clu, clustersToUpdate, fSkipRaw)
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
        nSamplesRaw = S0.traceDims(1);
        trWav_raw_clu = zeros(nSamplesRaw, nSitesSpike, nClusters, 'single');
        [tmrWav_raw_clu, tmrWav_raw_lo_clu, tmrWav_raw_hi_clu] = deal(zeros(nSamplesRaw, nSitesProbe, nClusters, 'single'));
    else
        [trWav_raw_clu, tmrWav_raw_clu, tmrWav_raw_lo_clu, tmrWav_raw_hi_clu] = deal([]);
    end

    if ~isempty(clustersToUpdate) % specific clusters we want to update
        updateCluster = ismember(1:nClusters, clustersToUpdate)';
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

    % Compute spkwav
    filteredTraces = getSpikeWaveforms(P, 0);
    for iCluster = 1:nClusters
        if updateCluster(iCluster)
            [mrWav_clu1, clusterSites1] = clu_wav_(S_clu, filteredTraces, iCluster, S0);
            if isempty(mrWav_clu1)
                continue;
            end
            
            [tmrWav_spk_clu(:, clusterSites1, iCluster), trWav_spk_clu(:, :, iCluster)] = ...
                deal(bit2uV_(mrWav_clu1, P));
        end
        if fVerbose, fprintf('.'); end
    end %clu

    % Compute spkraw
    if ~fSkipRaw
        rawTraces = []; % clear memory
        rawTraces = getSpikeWaveforms(P, 1);

        for iCluster = 1:nClusters
            if updateCluster(iCluster)
                [mrWav_clu1, clusterSites1, mrWav_lo_clu1, mrWav_hi_clu1] = clu_wav_(S_clu, rawTraces, iCluster, S0);
                if isempty(mrWav_clu1)
                    continue;
                end

                [tmrWav_raw_clu(:, clusterSites1, iCluster), trWav_raw_clu(:, :, iCluster)] = deal(meanSubtract(mrWav_clu1) * P.uV_per_bit);
                if isempty(mrWav_lo_clu1) || isempty(mrWav_hi_clu1)
                    tmrWav_raw_lo_clu(:, clusterSites1, iCluster) = zeros(nSamplesRaw, numel(clusterSites1));
                    tmrWav_raw_hi_clu(:, clusterSites1, iCluster) = zeros(nSamplesRaw, numel(clusterSites1));
                else
                    tmrWav_raw_lo_clu(:, clusterSites1, iCluster) = meanSubtract(mrWav_lo_clu1) * P.uV_per_bit;
                    tmrWav_raw_hi_clu(:, clusterSites1, iCluster) = meanSubtract(mrWav_hi_clu1) * P.uV_per_bit;
                end                
            end

            if fVerbose, fprintf('.'); end
        end % cluster
    end

    tmrWav_clu = tmrWav_spk_clu; %meanSubtract after or before?

    % measure waveforms
    [vrVmin_clu, viSite_min_clu] = min(permute(min(trWav_spk_clu),[2,3,1]),[],1);
    vrVmin_clu = abs(vrVmin_clu(:));
    viSite_min_clu = viSite_min_clu(:);

    S_clu = struct_add_(S_clu, vrVmin_clu, viSite_min_clu, ...
        trWav_spk_clu, tmrWav_spk_clu, trWav_raw_clu, tmrWav_raw_clu, tmrWav_clu, ...
        tmrWav_raw_lo_clu, tmrWav_raw_hi_clu);

    if fVerbose
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    end
end % function
