%--------------------------------------------------------------------------
function detect_(P)
    global tnWav_raw tnWav_spk trFet_spk;
    runtime_detect = tic;
    % Clear memory (S0 is cleared)
    set(0, 'UserData', []);
    [tnWav_raw, tnWav_spk, trFet_spk] = deal([]);

    S0 = file2spk_(P);

    if get_set_(P, 'fRamCache', 1)
        tnWav_raw = load_bin_(strrep(P.prmFile, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
        tnWav_spk = load_bin_(strrep(P.prmFile, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
    end
    trFet_spk = load_bin_(strrep(P.prmFile, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    S0.mrPos_spk = spk_pos_(S0, trFet_spk);

    % measure time
    S0.runtime_detect = toc(runtime_detect);
    fprintf('Detection took %0.1fs for %s\n', S0.runtime_detect, P.prmFile);

    set(0, 'UserData', S0);
    save0_(strrep(P.prmFile, '.prm', '_jrc.mat'));
    delete_(strrep(P.prmFile, '.prm', '_log.mat')); %delete log file when detecting
end %func
