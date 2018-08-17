%--------------------------------------------------------------------------
function importKiloSort(rezFile, sessionName)
    if nargin < 2 || ~isLegalFilename(sessionName)
        sessionName = 'imported-kilosort-session';
        fprintf('data will be saved to files beginning with ''%s''', sessionName);
    end

    if isstruct(rezFile)
        rez = rezFile;
    else
        try
            rez = load(rezFile, '-mat', 'rez');
        catch err
            if startsWith(err.message, 'Unable to read MAT-file')
                error('rezFile must be a MAT file');
            else
                rethrow(err);
            end
        end

        if isempty(rez)
            error('rez is empty');
        end
    end

    if isfield(rez, 'rez')
        rez = rez.rez;
    end

    if isfield(rez, 'Wphy')
        rez.W = rez.Wphy;
    end

    rez.W = gather(rez.W);
    rez.U = gather(rez.U);
    rez.mu = gather(rez.mu);

    % import spike location, time and cluster assignment
    spikeTimes = rez.st3(:, 1);
    spikeTemplates = rez.st3(:, 2);

    if size(rez.st3, 2) == 5
        spikeClusters = 1 + rez.st3(:, 5); % post-merging result
    else
        spikeClusters = spikeTemplates; % template/cluster
    end

    clusters = unique(spikeClusters);
    
    nTemplates = size(rez.simScore, 1);
    nClusters = numel(clusters);
    simScore = zeros(nClusters);
    
    clusterTemplates = cell(nClusters, 1);
    for iCluster = 1:nClusters
        cluster = clusters(iCluster);
        spikeClusterIndices = (spikeClusters == cluster); % spike indices for this cluster
        iClusterTemplates = unique(spikeTemplates(spikeClusterIndices)); % unique template indices for spikes in this cluster
        clusterTemplates{iCluster} = iClusterTemplates;
    end

    % consolidate cluster assignments
    for iCluster = 1:nClusters
        cluster = clusters(iCluster);
        spikeClusterIndices = (spikeClusters == cluster); % spike indices for this cluster
        
        if clusters(iCluster) ~= iCluster
            spikeClusters(spikeClusterIndices) = iCluster;
        end
        
        iClusterTemplates = clusterTemplates{iCluster}; % template IDs for this cluster

        % compute cluster sim score, Phy style
        sims = max(rez.simScore(iClusterTemplates, :), [], 1);
        
        for jCluster=iCluster:nClusters
            jClusterTemplates = clusterTemplates{jCluster};
%             if clusters(jCluster) <= nTemplates
%                 simScore(iCluster, jCluster) = sims(clusters(jCluster));
%             else
                simScore(iCluster, jCluster) = max(sims(jClusterTemplates));
%             end
            simScore(jCluster, iCluster) = simScore(iCluster, jCluster);
        end
    end
    
    % save the old clusters with gaps in them
    clustersGapped = clusters;
    clusters = (1:nClusters)';

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

    % compute the weighted average template for a given cluster and pick
    % its min site
    clusterSites = zeros('like', clusters);
    for iCluster = 1:nClusters
        iClusterTemplates = clusterTemplates{iCluster}; % templates for this cluster

        avgTemplate = squeeze(templates(iClusterTemplates, :, :));
        if numel(iClusterTemplates) > 1
            freqs = histcounts(spikeTemplates(spikeClusters == iCluster), numel(iClusterTemplates));
            weights = freqs/sum(freqs);
            t = zeros(size(avgTemplate, 2), size(avgTemplate, 3));
            for iWeight = numel(weights)
                t = t + weights(iWeight)*squeeze(avgTemplate(iWeight, :, :));
            end
            
            avgTemplate = t;
        end
        
        sampleMin = min(avgTemplate, [], 1);
        [~, clusterSites(iCluster)] = min(sampleMin); % cluster location
    end
    spikeSites = clusterSites(spikeClusters);

    % construct P from scratch
    P = struct();
    P.paramFile = [sessionName '.prm'];
    P.vcFile = rez.ops.fbinary;
    P.sampleRateHz = rez.ops.fs;

    % probe
    P.chanMap = rez.ops.chanMap;
    P.mrSiteXY = [rez.xcoords(:) rez.ycoords(:)];
    P.nChans = rez.ops.NchanTOT;
    P.viShank_site = rez.ops.kcoords(:)';
    P.vrSiteHW = [12 12]; % TODO: address this
    if isfield(rez, 'connected')
        P.chanMap = P.chanMap(rez.connected);
        P.mrSiteXY = P.mrSiteXY(rez.connected(:), :);
    end
    P.viChan_aux = setdiff(1:P.nChans, 1:max(P.chanMap));

    % Phy filter settings
    P.vcFilter = 'bandpass';
    P.freqLim = [500, .475 * P.sampleRateHz];

    P.blank_thresh = 0;
    P.dataType = 'int16'; % KS default
    P.fImportKilosort = 1;
    P.fTranspose_bin = 1;
    P.feature = 'kilosort';
    P.maxSite = 6.5; % TODO: allow user to set
    P.nSites_ref = 0; % TODO: address

    P.spkLim = ceil(size(templates, 2)/2) * [-1, 1];
    P.spkLim_raw = P.spkLim; % TODO: address

    P.corrLim = [.75 1];
    P.fDrift_merge = 0; % do not attempt drift correction
    P.nPaddingSamples = 100; % default
    P.qqFactor = 5; % default
    P.nPcPerChan = 1; % default
    P.nTime_clu = 1; % spikes detected and clustered over entire time series
    P.uV_per_bit = 1; % set this to unit scaling and deal with later
    P.spkRefrac_ms = .25; % default

    P.viSiteZero = [];
    P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site);
    P.useGPU = double(gpuDeviceCount() > 0);

%     S0 = kilosort2jrc_(P, int32(spikeTimes), int32(spikeSites));
    S0 = file2spk_(P, int32(spikeTimes), int32(spikeSites));
    P = saveProbe([sessionName '-probe.mat'], P);
    S0.P = P;

    % extract features
    spikeFeatures = permute(rez.cProjPC, [3 2 1]);
    fidFeatures = fopen([sessionName '_features.bin'], 'w');
    fwrite_(fidFeatures, spikeFeatures);
    fclose(fidFeatures);

    featureDims = size(spikeFeatures);
    S0.featureDims = featureDims;
    S0.rez = rez;
    
    set(0, 'UserData', S0);

    % construct S_clu from scratch
    S_clu = struct();
    % TODO: handle numel(clusters) ~= max(clusters);
    S_clu.nClusters = nClusters;
    S_clu.spikeClusters = spikeClusters;
    S_clu.clustersGapped = clustersGapped;
    S_clu.spikeClustersAuto = spikeTemplates;
    S_clu.clusterNotes = cell(nClusters, 1);
    S_clu.clusterSites = clusterSites;
    S_clu.simScore = simScore;

    S_clu.spikesByCluster = cell(1, nClusters);
    for iCluster = 1:nClusters
        S_clu.spikesByCluster{iCluster} = find(spikeClusters == iCluster);
    end

%     S0.mrPos_spk = spk_pos_(S0, spikeFeatures);
%     set(0, 'UserData', S0);
%
%     % cluster and describe
%     S0.S_clu = cluster_spacetime_(S0, P);
%     S0.S_clu = S_clu_new_(spikeClusters, S0);
%     S0.S_clu = S_clu_sort_(S0.S_clu, 'clusterSites');

    S0.S_clu = S_clu;

    % Save
    set(0, 'UserData', S0);

    exportParams(P.paramFile, [], 0);
    save0_([sessionName '_jrc.mat']);
end %func
