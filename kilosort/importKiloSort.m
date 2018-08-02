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

    global tnWav_raw tnWav_spk trFet_spk
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing
    spikeTimes = rez.st3(:, 1);

    if size(rez.st3, 2) == 5
        spikeClusters = 1 + rez.st3(:, 5); % post-merging result
    else
        spikeClusters = rez.st3(:, 2); % template/cluster
    end

    if isfield(rez, 'Wraw')
        Wraw = gather(rez.Wraw); % nC, nT, nClu
    else % reconstruct Wraw
        Wraw = zeros(size(rez.U, 1), size(rez.W, 1), numel(rez.mu));

        Urot = gather(rez.U);
        W = gather(rez.W);
        WrotInv = gather(inv(rez.Wrot));
        mu = gather(rez.mu);

        for k = 1:size(Urot, 3)
            Urot(:, :, k)  = WrotInv' * Urot(:,:,k);
        end

        for n = 1:size(rez.U, 2)
            Wraw(:,:,n) = mu(n) * squeeze(Urot(:,n,:)) * squeeze(W(:,n,:))';
        end
    end

    Wraw = -abs(Wraw);
    Wraw = permute(Wraw, [2,1,3]);
    mnMin_clu = squeeze_(min(Wraw, [], 1));
    [~, clusterSites] = min(mnMin_clu, [], 1); % cluster location

    % construct P from scratch
    P = struct();
    P.chanMap = rez.ops.chanMap;

    clusterSites = P.chanMap(clusterSites);
    viSite_spk = clusterSites(spikeClusters);

    return; % WIP

    S0 = file2spk_(P, int32(spikeTimes), int32(viSite_spk));
    S0.P = P;
    S0.S_ksort = S_ksort;
    tnWav_raw = load_bin_(strrep(P.rezFile, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
    tnWav_spk = load_bin_(strrep(P.rezFile, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
    trFet_spk = load_bin_(strrep(P.rezFile, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    S0.mrPos_spk = spk_pos_(S0, trFet_spk);
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
