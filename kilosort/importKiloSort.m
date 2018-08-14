%--------------------------------------------------------------------------
function importKiloSort(rezFile)
    % importKiloSort(rezFile)
    if isstruct(rezFile)
        rez = rezFile;
    else
        try
            rez = load(rezFile, '-mat', 'rez');
        catch err
            if startsWith(err.message, 'Unable to read MAT-file')
                error('rezFile must be a MAT file');
            else
                throw err;
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
    templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels

    sampleMin = squeeze(min(templates, [], 2));
    [~, clusterSites] = min(sampleMin, [], 2); % cluster location

    % construct P from scratch
    P = struct();
    P.paramFile = 'imported-kilosort-session.prm'; % TODO: allow user to set
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
%     P.nSites_ref = []; % JRC default; TODO: decide or allow user to set

    P.spkLim = ceil(size(templates, 2)/2) * [-1, 1];
    P.spkLim_raw = P.spkLim; % TODO: address

    P.viSiteZero = [];
    P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site);

    spikeSites = clusterSites(spikeTemplates);

    S0 = kilosort2jrc_(P, int32(spikeTimes), int32(spikeSites));
    S0.P = P;
    
    % extract features
    spikeFeatures = permute(rez.cProjPC, [2 3 1]);
    fidFeatures = fopen(strrep(P.paramFile, '.prm', '_features.bin'), 'W');
    fwrite_(fidFeatures, spikeFeatures);
    fclose(fidFeatures);
    
    featureDims = size(spikeFeatures);
    S0.featureDims = featureDims;
    S0.rez = rez;

%     S0.mrPos_spk = spk_pos_(S0, spikeFeatures);
%     set(0, 'UserData', S0);
% 
%     % cluster and describe
%     S0.S_clu = cluster_spacetime_(S0, P);
%     S0.S_clu = S_clu_new_(spikeClusters, S0);
%     S0.S_clu = S_clu_sort_(S0.S_clu, 'clusterSites');
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.paramFile, '.prm', '_jrc.mat'));
end %func
