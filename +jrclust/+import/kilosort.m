function hClust = kilosort(rezFile)
    %KILOSORT Import a Kilosort session from rez.mat
    rezFile_ = jrclust.utils.absPath(rezFile);
    if isempty(rezFile_)
        error('Could not find file ''%s''', rezFile);
    end

    rezFile = rezFile_;

    try
        load(rezFile, 'rez');
    catch ME
        error('Failed to load ''%s'': %s', rezFile, ME.message);
    end

    rezFields = fieldnames(rez); %#ok<NODEF>
    for i = 1:numel(rezFields)
        fn = rezFields{i};
        if isa(rez.(fn), 'gpuArray')
            rez.(fn) = jrclust.utils.tryGather(rez.(fn));
        end
    end

    spikeTimes = rez.st3(:, 1);
    spikeTemplates = rez.st3(:, 2);

    amplitudes = rez.st3(:, 3);

    if size(rez.st3, 2) > 4
        spikeClusters = 1 + rez.st3(:, 5);
    else
        spikeClusters = spikeTemplates;
    end

    [clusterIDs, ~, indices] = unique(spikeClusters);
    goodClusters = clusterIDs(clusterIDs > 0);
    junkClusters = setdiff(clusterIDs, goodClusters);
    clusterIDsNew = [junkClusters' 1:numel(goodClusters)]';
    spikeClusters = clusterIDsNew(indices);
    nTemplates = size(rez.simScore, 1);

    nClusters = numel(goodClusters);

    clusterTemplates = arrayfun(@(iCluster) unique(spikeTemplates(spikeClusters == iCluster)), 1:nClusters, 'UniformOutput', 0);

    % compute templates
    nt0 = size(rez.W, 1);
    U = rez.U;
    W = rez.W;

    Nfilt = size(W, 2);
    Nchan = rez.ops.Nchan;

    templates = zeros(Nchan, nt0, Nfilt, 'single');
    for iNN = 1:size(templates,3)
       templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
    end
    templates = -abs(permute(templates, [3 2 1])); % nTemplates x nSamples x nChannels

    % compute the weighted average template for a given cluster and pick its min site
    clusterSites = zeros('like', clusterIDs);
    for iCluster = 1:nClusters
        iTemplates = clusterTemplates{iCluster}; % templates for this cluster

        meanTemplate = squeeze(templates(iTemplates, :, :));
        if numel(iTemplates) > 1
            freqs = histcounts(spikeTemplates(spikeClusters == iCluster), numel(iTemplates));
            weights = freqs/sum(freqs); % weight sum of templates by frequency of occurrence in this cluster
            t = zeros(size(meanTemplate, 2), size(meanTemplate, 3));

            for iWeight = numel(weights)
                t = t + weights(iWeight)*squeeze(meanTemplate(iWeight, :, :));
            end

            meanTemplate = t;
        end

        sampleMin = min(meanTemplate, [], 1);
        [~, clusterSites(iCluster)] = min(sampleMin); % cluster location
    end
    spikeSites = clusterSites(spikeClusters);

    dRes = struct('spikeTimes', spikeTimes, ...
                  'spikeSites', spikeSites);
    sRes = struct('spikeClusters', spikeClusters, ...
                  'spikeTemplates', spikeTemplates);

    hCfg = jrclust.Config();
    
end

