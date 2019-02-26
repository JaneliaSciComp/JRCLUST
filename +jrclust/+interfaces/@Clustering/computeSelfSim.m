function selfSim = computeSelfSim(obj, iCluster)
    %COMPUTESELFSIM Get similarity between bottom and top half (vpp-wise) of cluster
    if nargin < 2
        iCluster = [];
    end

    if isempty(iCluster)
        obj.hCfg.updateLog('selfSim', 'Computing cluster self-similarity', 1, 0);

        selfSim = zeros(1, obj.nClusters);
        for iCluster = 1:obj.nClusters
            selfSim(iCluster) = iSelfSim(obj, iCluster);
        end

        obj.hCfg.updateLog('selfSim', 'Finished computing cluster self-similarity', 0, 1);
    else
        selfSim = iSelfSim(obj, iCluster);
    end
end

%% LOCAL FUNCTIONS
function selfSim = iSelfSim(hClust, iCluster)
    %ISELFSIM Kernel of the self-similarity method (for a given cluster)
    %   low means bad. return 1-corr score
    MAX_SAMPLE = 4000;

    clusterSpikes_ = hClust.getCenteredSpikes(iCluster);
    clusterSpikes_ = jrclust.utils.subsample(clusterSpikes_, MAX_SAMPLE);

    centeredWf = hClust.spikesRaw(:, :, clusterSpikes_);
    centeredVpp = squeeze(squeeze(max(centeredWf(:, 1, :)) - min(centeredWf(:, 1, :))));
    [~, argsort] = sort(centeredVpp);

    imid = round(numel(argsort)/2);

    % compute correlation between low-vpp waveforms and high-vpp waveforms
    % only on spikes occurring on the center site
    lowHalf = jrclust.utils.meanSubtract(mean(centeredWf(:, :, argsort(1:imid)), 3));
    highHalf = jrclust.utils.meanSubtract(mean(centeredWf(:, :, argsort(imid+1:end)), 3));

    selfSim = colZScore(lowHalf(:), highHalf(:));
end

function C = colZScore(X, Y)
    % https://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab
    % compute Z-scores for columns of X and Y
    X = bsxfun(@minus, X, mean(X)); % zero-mean
    X = bsxfun(@times, X, 1./sqrt(sum(X.^2))); %% L2-normalization

    Y = bsxfun(@minus, Y, mean(Y)); % zero-mean
    Y = bsxfun(@times, Y, 1./sqrt(sum(Y.^2))); %% L2-normalization

    C = X' * Y;
end
