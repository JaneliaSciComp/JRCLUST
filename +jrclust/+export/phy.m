function phy(hClust)
    %PHY Export JRCLUST data to .npy format for Phy
    if exist('writeNPY', 'file') ~= 2
        warning('Please make sure you have npy-matlab installed and on your path (https://github.com/kwikteam/npy-matlab)');
        return;
    elseif ~isa(hClust, 'jrclust.sort.DensityPeakClustering')
        error('Phy export not supported for this type of clustering.');
    end

    hCfg = hClust.hCfg;

    [nSites, ~, nSpikes] = size(hClust.spikeFeatures);

    % spikeAmps/amplitudes
    ampFile = fullfile(hCfg.outputDir, 'amplitudes.npy');
    writeNPY(abs(hClust.spikeAmps), ampFile);
    hCfg.updateLog('phy', sprintf('Saved spikeAmps to %s', ampFile), 0, 0);

    % spikeTimes/spike_times
    timeFile = fullfile(hCfg.outputDir, 'spike_times.npy');
    writeNPY(uint64(hClust.spikeTimes), timeFile);
    hCfg.updateLog('phy', sprintf('Saved spikeTimes to %s', timeFile), 0, 0);

    % spikeSites/spike_sites (0 offset)
    siteFile = fullfile(hCfg.outputDir, 'spike_sites.npy');
    writeNPY(int32(hClust.spikeSites) - 1, siteFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved spikeSites to %s', siteFile), 0, 0);

    % siteMap/channel_map (0 offset)
    mapFile = fullfile(hCfg.outputDir, 'channel_map.npy');
    writeNPY((int32(hCfg.siteMap) - 1)', mapFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved siteMap to %s', mapFile), 0, 0);

    % siteLoc/channel_positions
    locFile = fullfile(hCfg.outputDir, 'channel_positions.npy');
    writeNPY(hCfg.siteLoc, locFile);
    hCfg.updateLog('phy', sprintf('Saved siteLoc to %s', locFile), 0, 0);

    % similar templates
    simFile = fullfile(hCfg.outputDir, 'similar_templates.npy');
    writeNPY(single(hClust.waveformSim), simFile);
    hCfg.updateLog('phy', sprintf('Saved waveformSim to %s', simFile), 0, 0);

    % template feature ind
    tfiFile = fullfile(hCfg.outputDir, 'template_feature_ind.npy');
    [~, argsort] = sort(hClust.waveformSim, 1);
    tfi = argsort((hClust.nClusters-5):end, 1:end)';
    writeNPY(uint32(tfi), tfiFile);
    hCfg.updateLog('phy', sprintf('Saved template feature inds to %s', tfiFile), 0, 0);

    % template features
    tfFile = fullfile(hCfg.outputDir, 'template_features.npy');

    % meanWfLocalRaw, spikesRaw
    tf = zeros(hClust.nSpikes, 6, 'single');
    for i = 1:hClust.nClusters
        near = tfi(i, :);
        iSpikes = hClust.spikesByCluster{i};
        waves = reshape(hClust.spikesFilt(:, :, iSpikes), size(hClust.spikesFilt, 1)*size(hClust.spikesFilt, 2), numel(iSpikes));

        for j = 1:numel(near)
            wave_mu = reshape(hClust.meanWfLocal(:, :, j), 1, size(hClust.meanWfLocal, 1)*size(hClust.meanWfLocal,2));
            template_features = (single(wave_mu)*single(waves))./(numel(wave_mu).*std(single(wave_mu)).*std(single(waves), [], 1));
            tf(iSpikes,j)=template_features;
        end
    end
    writeNPY(tf, tfFile);
    hCfg.updateLog('phy', sprintf('Saved template features to %s', tfFile), 0, 0);

    % templates (mwf)
    templatesFile = fullfile(hCfg.outputDir, 'templates.npy');
    templates = zeros(hClust.nClusters, diff(hClust.hCfg.evtWindowSamp) + 1, numel(hClust.siteRMS), 'single');
    for i = 1:hClust.nClusters
        templates(i, :, :) = hClust.meanWfGlobal(:, :, i)./max(max(abs(hClust.meanWfGlobal(:, :, i))))./max(max(abs(hClust.meanWfGlobal(:,:,i))));
    end
    writeNPY(templates, templatesFile);
    hCfg.updateLog('phy', sprintf('Saved templates to %s', templatesFile), 0, 0);
    
    % templates_ind (meaningless?)
    tiFile = fullfile(hCfg.outputDir, 'templates_ind.npy');
    nsites = numel(hClust.siteRMS);
    templatesInd = zeros(hClust.nClusters,nsites);
    for i = 1:hClust.nClusters
        templatesInd(i, :, :) = (1:nsites) - 1;
    end
    writeNPY(templatesInd, tiFile);
    hCfg.updateLog('phy', sprintf('Saved templates to %s', tiFile), 0, 0);
           
    % spike_templates
    stFile = fullfile(hCfg.outputDir, 'spike_templates.npy');
    st = hClust.spikeClusters - 1;
    st(st < 0) = 0;
    writeNPY(uint32(st), stFile);
    hCfg.updateLog('phy', sprintf('Saved templates to %s', stFile), 0, 0);
        
    % whitening mat
    wmFile = fullfile(hCfg.outputDir, 'whitening_mat.npy');
    nsites = numel(hClust.siteRMS);
    writeNPY(eye(nsites, nsites).*.4, wmFile);
    hCfg.updateLog('phy', sprintf('Saved whitening mat to %s', wmFile), 0, 0);

    % whitening mat inv
    wmiFile = fullfile(hCfg.outputDir, 'whitening_mat_inv.npy');
    nsites = numel(hClust.siteRMS);
    writeNPY(eye(nsites, nsites).*.4, wmiFile);
    hCfg.updateLog('phy', sprintf('Saved whitening mat inv to %s', wmiFile), 0, 0);
    
    % spikeClusters/spike_clusters
    spikeClusters = hClust.spikeClusters;
    spikeClusters(spikeClusters < 0) = 0;
    clusterFile = fullfile(hCfg.outputDir, 'spike_clusters.npy');
    writeNPY(uint32(spikeClusters) - 1, clusterFile); % -1 for zero indexing
    hCfg.updateLog('phy', sprintf('Saved spikeClusters to %s', clusterFile), 0, 0);

    % write out PC features
    if ismember(hCfg.clusterFeature, {'pca', 'gpca'})
        % take just primary peak
        spikeFeatures = squeeze(hClust.spikeFeatures(:, 1, :))'; % nSpikes x nFeatures
        if hCfg.nPCsPerSite == 1
            pcFeatures = reshape(spikeFeatures, size(spikeFeatures, 1), 1, []);
            pcFeatures(:,2,:) = 0;
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
        pfi = zeros(hClust.nClusters, nSites, 'uint32');
        for i = 1:hClust.nClusters
            pfi(i, :) = hCfg.siteNeighbors(1:nSites, hClust.clusterSites(i)) - 1; % - 1 for zero indexing
        end
        writeNPY(pfi, indFile);
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