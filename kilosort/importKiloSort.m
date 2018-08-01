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
    end

    if isempty(rez)
        error('rez is empty');
    end

    global tnWav_raw tnWav_spk trFet_spk
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing

    % construct P from scratch
    spikeTimes = S_ksort.rez.st3(:,1); %spike time
    if get_set_(P, 'fMerge_post', 0) && size(S_ksort.rez.st3, 2) == 5
        viClu = 1 + S_ksort.rez.st3(:,5); % post-merging result
    else
        viClu = S_ksort.rez.st3(:,2); % template/cluster
    end

    tnWav_clu = S_ksort.rez.Wraw; %nC, nT, nClu
    tnWav_clu = -abs(tnWav_clu);
    tnWav_clu = permute(tnWav_clu, [2,1,3]);
    mnMin_clu = squeeze_(min(tnWav_clu, [], 1));
    [~, viSite_clu] = min(mnMin_clu, [], 1); %cluster location
    viSite = 1:numel(P.chanMap);
    viSite(P.viSiteZero) = [];
    viSite_clu = viSite(viSite_clu);
    viSite_spk = viSite_clu(viClu);

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
    S0.S_clu = S_clu_new_(viClu, S0);
    S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu');
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.rezFile, '.prm', '_jrc.mat'));
    describe_(S0);
end %func
