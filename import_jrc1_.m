%--------------------------------------------------------------------------
function import_jrc1_(vcFile_prm)
    % import jrc1 so that i can visualize the output
    global tnWav_raw tnWav_spk trFet_spk
    % convert jrc1 format (_clu and _evt) to jrc3 format. no overwriting
    % receive spike location, time and cluster number. the rest should be taken care by jrc3 processing

    % Load info from previous version: time, site, spike
    P = loadParam_(vcFile_prm);
    Sevt = load(strrep(P.vcFile_prm, '.prm', '_evt.mat'));
    if isfield(Sevt, 'Sevt'), Sevt = Sevt.Sevt; end
    S0 = struct('viTime_spk', Sevt.viSpk, 'viSite_spk', Sevt.viSite, 'P', P);

    [tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, S0.viTime_spk, S0.viSite_spk);
    set(0, 'UserData', S0);

    % Save to file
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk);
    write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk);
    save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));

    % cluster and describe
    sort_(P);
    S0 = get(0, 'UserData');
    Sclu = load_(strrep(P.vcFile_prm, '.prm', '_clu.mat'));
    if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end
    if ~isempty(Sclu)
        S0.S_clu.viClu = Sclu.viClu; %skip FigRD step for imported cluster
    end
    set(0, 'UserData', S0);

    describe_(P.vcFile_prm);
end %func
