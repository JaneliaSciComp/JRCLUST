function [retained, spikeFeatures] = doAutoSplit(sampledSpikes, spikeFeatures, nSplits, hCfg)
    %DOAUTOSPLIT
    % TODO: ask users number of clusters and split multi-way
    %Make automatic split of clusters using PCA + kmeans clustering
    [~, pcaFeatures, ~] = pca(double(sampledSpikes'), 'Centered', 1, 'NumComponents', 3);
    spikeFeatures = double([spikeFeatures pcaFeatures]);

    if hCfg.getOr('fUseMikeSplit', 0)
        inClust = MikeSplit(sampledSpikes, spikeFeatures, nSplits);
    end

    spikeFeatures = pcaFeatures;
    nSpikes = size(sampledSpikes, 2);

    % ask how many clusters there are
    try
        % kmean clustering into 2
        kmAssigns = kmeans(spikeFeatures, nSplits);
        d12 = madDist(spikeFeatures(kmAssigns == 1, :)', spikeFeatures(kmAssigns == 2, :)');

        if hCfg.verbose
            fprintf('MAD distance: %f\n', d12);
        end

        if hCfg.getOr('fUseMikeSplit', 0)
            retained = inClust;
        else
            retained = (kmAssigns == 1);
        end
    catch % not enough features for k-means to automatically split
        retained = false(nSpikes, 1);
        retained(1:end/2) = 1;
    end
end

%% LOCAL FUNCTIONS
function d12 = madDist(features1, features2)
    % distance between two clusters
    if ~ismatrix(features1)
        features1 = reshape(features1, [], size(features1, 3));
    end
    if ~ismatrix(features2)
        features2 = reshape(features2, [], size(features2, 3));
    end

    med1 = median(features1, 2);
    med2 = median(features2, 2);

    medDiffs = med1 - med2;
    norm12 = norm(medDiffs, 2);
    vrFet12_med1 = medDiffs/norm12;

    mad1 = median(abs(vrFet12_med1' * bsxfun(@minus, features1, med1)));
    mad2 = median(abs(vrFet12_med1' * bsxfun(@minus, features2, med2)));

    d12 = norm12 / norm([mad1, mad2], 2);
end