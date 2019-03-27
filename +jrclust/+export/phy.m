function phy(hCfg)
    %PHY Export JRCLUST data to .npy format for Phy
    if exist('writeNPY', 'file') ~= 2
        warning('Please make sure you have npy-matlab installed (https://github.com/kwikteam/npy-matlab)');
        return;
    end

    res = load(hCfg.resFile);

    if ~isfield(res, 'spikeTimes')
        return;
    end

    % spikeAmps/amplitudes
    ampFile = fullfile(hCfg.outputDir, 'amplitudes.npy');
    writeNPY(uint64(res.spikeAmps), ampFile);
    hCfg.updateLog('phy', sprintf('Saved spikeAmps to %s', ampFile), 0, 0);

    % spikeTimes/spike_times
    timeFile = fullfile(hCfg.outputDir, 'spike_times.npy');
    writeNPY(uint64(res.spikeTimes), timeFile);
    hCfg.updateLog('phy', sprintf('Saved spikeTimes to %s', timeFile), 0, 0);

    % spikeSites/spike_sites (0 offset)
    siteFile = fullfile(hCfg.outputDir, 'spike_sites.npy');
    writeNPY(int32(res.spikeSites) - 1, siteFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved spikeSites to %s', siteFile), 0, 0);

    % siteMap/channel_map (0 offset)
    mapFile = fullfile(hCfg.outputDir, 'channel_map.npy');
    writeNPY(int32(hCfg.siteMap) - 1, mapFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved siteMap to %s', mapFile), 0, 0);

    % siteLoc/channel_positions
    locFile = fullfile(hCfg.outputDir, 'channel_positions.npy');
    writeNPY(hCfg.siteLoc, locFile);
    hCfg.updateLog('phy', sprintf('Saved siteLoc to %s', locFile), 0, 0);

    % spikeClusters/spike_clusters
    if isfield(res, 'hClust') && ~isempty(res.hClust.spikeClusters)
        spikeClusters = res.hClust.spikeClusters;
    elseif isfield(res, 'spikeClusters')
        spikeClusters = res.spikeClusters;
    else
        spikeClusters = [];
    end
    spikeClusters(spikeClusters < 0) = 0;

    clusterFile = fullfile(hCfg.outputDir, 'spike_clusters.npy');
    writeNPY(uint32(spikeClusters) - 1, clusterFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved spikeClusters to %s', clusterFile), 0, 0);

    % write out PC features
    if ismember(hCfg.clusterFeature, {'pca', 'gpca'})
        fid = fopen(hCfg.featuresFile, 'r');
        spikeFeatures = reshape(fread(fid, inf, '*single'), res.featuresShape);
        fclose(fid);

        [nSites, nSpikes] = deal(res.featuresShape(1), res.featuresShape(3));

        % take just primary peak
        spikeFeatures = squeeze(spikeFeatures(:, 1, :))'; % nSpikes x nFeatures
        if hCfg.nPCsPerSite == 1
            pcFeatures = reshape(spikeFeatures, size(spikeFeatures, 1), 1, []);
        else
            nSites = nSites/hCfg.nPCsPerSite;
            pcFeatures = zeros(nSpikes, hCfg.nPCsPerSite, nSites, 'single');
            for i = 1:hCfg.nPCsPerSite
                pcFeatures(:, i, :) = spikeFeatures(:, ((i-1)*nSites+1):i*nSites);
            end
        end

        featuresFile = fullfile(hCfg.outputDir, 'pc_features.npy');
        writeNPY(pcFeatures, featuresFile);
        hCfg.updateLog('phy', sprintf('Saved spikeFeatures to %s', featuresFile), 0, 0);

        indFile = fullfile(hCfg.outputDir, 'pc_feature_ind.npy');
        writeNPY(uint32(hCfg.siteNeighbors(1:nSites, res.spikeSites)') - 1, indFile); % -1 for zero indexing
        hCfg.updateLog('phy', sprintf('Saved spikeFeature indices to %s', indFile), 0, 0);
    end

    % param file
    if exist(fullfile(hCfg.outputDir, 'params.py'), 'file') ~= 2
        paramFile = fullfile(hCfg.outputDir, 'params.py');
        fid = fopen(paramFile, 'w');

        rawRecordings = cellfun(@(x) sprintf('r''%s''', x), hCfg.rawRecordings, 'UniformOutput', 0);
        rawRecordings = ['[' strjoin(rawRecordings, ', ') ']'];
        fprintf(fid, 'dat_path = %s\n', rawRecordings);

        fprintf(fid, 'n_channels_dat = %i\n', hCfg.nChans);

        fprintf(fid, 'dtype = ''%s''\n', dtype2NPY(hCfg.dataType));

        fprintf(fid, 'offset = %d\n', hCfg.headerOffset);

        if hCfg.sampleRate ~= floor(hCfg.sampleRate)
            fprintf(fid,'sample_rate = %i\n', hCfg.sampleRate);
        else
            fprintf(fid,'sample_rate = %i.\n', hCfg.sampleRate);
        end

        fprintf(fid,'hp_filtered = False');
        fclose(fid);

        hCfg.updateLog('phy', sprintf('Saved params to %s', paramFile), 0, 0);
    end
end

%% LOCAL FUNCTIONS
function dtype = dtype2NPY(dtype)
    switch dtype
        case 'single'
            dtype = 'float32';
        case 'double'
            dtype = 'float64';
    end
end
