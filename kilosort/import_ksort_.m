%--------------------------------------------------------------------------
function import_ksort_(vcFile_prm)
    % import_ksort_(P, fSort)
    % import_ksort_(vcFile_prm, fSort)
    % fMerge_post = 0;
    % import kilosort result
    if isstruct(vcFile_prm)
        P = vcFile_prm;
        vcFile_prm = P.prmFile;
    else
        P = loadParams(vcFile_prm); %makeParam_kilosort_
    end

    % S_ksort = load(strrep(P.prmFile, '.prm', '_ksort.mat')); % contains rez structure
    global tnWav_raw tnWav_spk trFet_spk
    % convert jrc1 format (_clu and _evt) to jrc format. no overwriting
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing

    % Create a prm file to start with. set the filter parameter correctly. features?
    if isempty(P), return; end
    try
        S_ksort = load(strrep(P.prmFile, '.prm', '_ksort.mat')); %get site # and
    catch % import rez.mat -- acl
        rez = load('rez.mat', 'rez');
        S_ksort.rez = rez.rez; % -_-
        S_ksort.P = P;
        S_ksort.runtime_ksort = 0; % don't have this
        struct_save_(S_ksort, strrep(vcFile_prm, '.prm', '_ksort.mat'), 1);
    end
    viTime_spk = S_ksort.rez.st3(:,1); %spike time
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
    viSite = 1:numel(P.viSite2Chan);
    viSite(P.viSiteZero) = [];
    viSite_clu = viSite(viSite_clu);
    viSite_spk = viSite_clu(viClu);

    S0 = file2spk_(P, int32(viTime_spk), int32(viSite_spk));
    S0.P = P;
    S0.S_ksort = S_ksort;
    tnWav_raw = load_bin_(strrep(P.prmFile, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
    tnWav_spk = load_bin_(strrep(P.prmFile, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
    trFet_spk = load_bin_(strrep(P.prmFile, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    S0.mrPos_spk = spk_pos_(S0, trFet_spk);
    set(0, 'UserData', S0);

    % cluster and describe
    S0.S_clu = cluster_spacetime_(S0, P);
    S0.S_clu = S_clu_new_(viClu, S0);
    S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu');
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.prmFile, '.prm', '_jrc.mat'));
    describe_(S0);
end %func
