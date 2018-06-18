%--------------------------------------------------------------------------
function import_ksort_(vcFile_prm, fSort)
    % import_ksort_(P, fSort)
    % import_ksort_(vcFile_prm, fSort)
    % fMerge_post = 0;
    % import kilosort result
    % fSort: do sort using jrc
    if isstruct(vcFile_prm)
        P = vcFile_prm;
        vcFile_prm = P.vcFile_prm;
    else
        P = loadParam_(vcFile_prm); %makeParam_kilosort_
    end

    % S_ksort = load(strrep(P.vcFile_prm, '.prm', '_ksort.mat')); % contains rez structure
    if nargin<2, fSort = 0; end %
    global tnWav_raw tnWav_spk trFet_spk
    % convert jrc1 format (_clu and _evt) to jrc format. no overwriting
    % receive spike location, time and cluster number. the rest should be taken care by jrc processing

    % Create a prm file to start with. set the filter parameter correctly. features?
    if isempty(P), return; end
    S_ksort = load(strrep(P.vcFile_prm, '.prm', '_ksort.mat')); %get site # and
    viTime_spk = S_ksort.rez.st3(:,1) - 6; %spike time (apply shift factor)
    viClu = S_ksort.rez.st3(:,2); % cluster

    viClu_post = 1 + S_ksort.rez.st3(:,5); %post-merging result
    tnWav_clu = S_ksort.rez.Wraw; %nC, nT, nClu
    tnWav_clu = -abs(tnWav_clu);
    tnWav_clu = permute(tnWav_clu, [2,1,3]);
    mnMin_clu = squeeze_(min(tnWav_clu, [], 1));
    [~, viSite_clu] = min(mnMin_clu, [], 1); %cluster location
    viSite = 1:numel(P.viSite2Chan);
    viSite(P.viSiteZero) = [];
    viSite_clu = viSite(viSite_clu);
    viSite_spk = viSite_clu(viClu);
    % vnAmp_spk = S_ksort.rez.st3(:,3);

    % S0 = struct('viTime_spk', int32(viTime_spk), 'viSite_spk', int32(viSite_spk), 'P', P, 'S_ksort', S_ksort);
    S0 = file2spk_(P, int32(viTime_spk), int32(viSite_spk));
    S0.P = P;
    S0.S_ksort = S_ksort;
    tnWav_raw = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
    tnWav_spk = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
    trFet_spk = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    S0.mrPos_spk = spk_pos_(S0, trFet_spk);
    set(0, 'UserData', S0);

    % cluster and describe
    S0.S_clu = cluster_spacetime_(S0, P);
    
    if get_set_(P, 'fMerge_post', 0)
        S0.S_clu = S_clu_new_(viClu_post, S0);
    else
        S0.S_clu = S_clu_new_(viClu, S0);
    end
    S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu');
    set(0, 'UserData', S0);

    % Save
    save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
    describe_(S0);
end %func
