function export_jrc1_(vcFile_prm)
    % Export to version 1 format (_clu.mat) and (_evt.mat)
    % error('export_jrc1_: not implemented yet');
    % Load info from previous version: time, site, spike
    P = loadParam_(vcFile_prm);
    S0 = load(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));

    fprintf('Exporting to JRCLUST ver.1 format\n');
    % Build Sevt and save
    Sevt = S0;
    Sevt = rmfield_(Sevt, 'S_clu', 'cmrFet_site', 'cvrTime_site', 'cvrVpp_site', ...
    'dimm_raw', 'dimm_spk', 'miSites_fet', 'mrFet', 'runtime_detect', 'runtime_sort', ...
    'viSite_spk', 'viT_offset_file', 'viTime_spk', 'vrAmp_spk');
    Sevt.cvrSpk_site = vr2cell_(S0.vrAmp_spk, S0.cviSpk_site); %S0.cvrVpp_site;
    % Sevt.miSites_fet = S0.miSites_fet;
    [Sevt.mrPv, Sevt.mrWav_spk, Sevt.trFet] = deal([]);
    Sevt.viSite = S0.viSite_spk;
    Sevt.viSpk = S0.viTime_spk;
    Sevt.vrSpk = S0.vrAmp_spk;
    Sevt.vrThresh_uV = bit2uV_(S0.vrThresh_site, S0.P);
    Sevt.dimm_fet = [1, S0.dimm_fet(:)'];
    write_struct_(strrep(P.vcFile_prm, '.prm', '_evt.mat'), Sevt);
    fprintf('\tEvent struct (Sevt) exported.\n\t');
    assignWorkspace_(Sevt);

    % S_clu
    if isfield(S0, 'S_clu')
        Sclu = S0.S_clu;
        Sclu = rmfield_(Sclu, 'tmrWav_raw_clu', 'tmrWav_spk_clu', 'trWav_raw_clu', 'trWav_spk_clu', 'viSite_min_clu', 'vrVmin_clu');
        Sclu.trWav_dim = [];
        Sclu.viSite = S0.viSite_spk;
        Sclu.viTime = S0.viTime_spk;
        Sclu.vrSnr_Vmin_clu = S0.S_clu.vrSnr_clu;
        write_struct_(strrep(P.vcFile_prm, '.prm', '_clu.mat'), Sclu);
        fprintf('\tCluster struct (Sclu) exported.\n\t');
        assignWorkspace_(Sclu);
    end
end %func
