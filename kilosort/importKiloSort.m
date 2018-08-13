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

    global spikeTraces spikeWaveforms spikeFeatures
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing
    spikeTimes = rez.st3(:, 1);

    if size(rez.st3, 2) == 5
        spikeClusters = 1 + rez.st3(:, 5); % post-merging result
    else
        spikeClusters = rez.st3(:, 2); % template/cluster
    end

    if isfield(rez, 'Wraw')
        unwhitenedTemplates = gather(rez.Wraw); % nChannels, nSamples, nClusters
    else % reconstruct Wraw
        unwhitenedTemplates = zeros(size(rez.U, 1), size(rez.W, 1), numel(rez.mu));

        Urot = gather(rez.U);
        W = gather(rez.W);
        WrotInv = gather(inv(rez.Wrot));
        mu = gather(rez.mu);

        for k = 1:size(Urot, 3)
            Urot(:, :, k)  = WrotInv' * Urot(:, :, k);
        end

        for n = 1:size(rez.U, 2)
            unwhitenedTemplates(:, :, n) = mu(n) * squeeze(Urot(:, n, :)) * squeeze(W(:, n, :))';
        end
    end

    unwhitenedTemplates = -abs(unwhitenedTemplates);
    unwhitenedTemplates = permute(unwhitenedTemplates, [2, 1, 3]);
    sampleMin = squeeze(min(unwhitenedTemplates, [], 1));
    [~, clusterSites] = min(sampleMin, [], 1); % cluster location

    % construct P from scratch
    P = struct();
    P.algorithm = 'KiloSort';
    P.chanMap = rez.ops.chanMap;
    if isfield(rez, 'connected')
        P.chanMap = P.chanMap(rez.connected);
    end
    P.dataType = 'int16'; % KS default
    P.fTranspose_bin = 1;
    P.feature = 'kilosort';
    P.maxSite = rez.ops.Nchan; % maybe -- acl
    P.nChans = rez.ops.Nchan;
    P.nSites_ref = []; % JRC default; TODO: decide or allow user to set
    P.paramFile = 'imported-kilosort-session.prm';
    P.vcFile = rez.ops.fbinary;

    clusterSites = P.chanMap(clusterSites);
    spikeSites = clusterSites(spikeClusters);

    S0 = file2spk_(P, int32(spikeTimes), int32(spikeSites));
    S0.P = P;

    return; % WIP

    spikeTraces = load_bin_(strrep(P.rezFile, '.prm', '_traces.bin'), 'int16', S0.traceDims);
    spikeWaveforms = load_bin_(strrep(P.rezFile, '.prm', '_waveforms.bin'), 'int16', S0.waveformDims);
    spikeFeatures = load_bin_(strrep(P.rezFile, '.prm', '_features.bin'), 'single', S0.featureDims);

    S0.mrPos_spk = spk_pos_(S0, spikeFeatures);
    set(0, 'UserData', S0);

    % cluster and describe
    S0.S_clu = cluster_spacetime_(S0, P);
    S0.S_clu = S_clu_new_(spikeClusters, S0);
    S0.S_clu = S_clu_sort_(S0.S_clu, 'clusterSites');
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.rezFile, '.prm', '_jrc.mat'));
    describe_(S0);
end %func
