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
    clusterTemplates = zeros('like', clusters);

    % consolidate cluster assignments
    for iCluster = 1:numel(clusters)
        spikeClusterIndices = (spikeClusters == clusters(iCluster));

        if clusters(iCluster) ~= iCluster
            spikeClusters(spikeClusterIndices) = iCluster;
        end
        
        % get the mode of templates for this cluster
        spikeClusterTemplates = spikeTemplates(spikeClusterIndices);
        clusterTemplates(iCluster) = mode(spikeClusterTemplates);
    end
    
    % save the old clusters with gaps in them
    clustersGapped = clusters;
    clusters = 1:numel(clusters);

    nClusters = clusters(end);

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
    templates = permute(templates, [3 2 1]); % nTemplates x nSamples x nChannels

    sampleMin = squeeze(min(templates, [], 2));
    [~, clusterSites] = min(sampleMin, [], 2); % cluster location
    clusterSites = clusterSites(clusterTemplates);

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

    P.corrLim = [.9 1]; % default
    P.fDrift_merge = 0; % do not attempt drift correction
    P.nTime_clu = 1; % spikes detected and clustered over entire time series
    P.uV_per_bit = 1; % set this to unit scaling and deal with later
    P. spkRefrac_ms = .25; % default

    P.viSiteZero = [];
    P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site);

    spikeSites = clusterSites(spikeClusters);

    S0 = kilosort2jrc_(P, int32(spikeTimes), int32(spikeSites));
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
    S_clu.simScore = rez.simScore(clusterTemplates, clusterTemplates);
    S_clu.clusterSites = clusterSites;

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
