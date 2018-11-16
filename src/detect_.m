%--------------------------------------------------------------------------
function detect_(P)
    global tnWav_raw tnWav_spk trFet_spk;
    runtime_detect = tic;
    % Clear memory (S0 is cleared)
    set(0, 'UserData', []);
    [tnWav_raw, tnWav_spk, trFet_spk] = deal([]);

    S0 = file2spk_(P);

    if get_set_(P, 'fRamCache', 1)
        tnWav_raw = jrclust.utils.readBin(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
        tnWav_spk = jrclust.utils.readBin(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
    end
    trFet_spk = jrclust.utils.readBin(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    S0.mrPos_spk = jrclust.utils.spikePos(S0.viSite_spk, trFet_spk, P);

    % measure time
    S0.runtime_detect = toc(runtime_detect);
    fprintf('Detection took %0.1fs for %s\n', S0.runtime_detect, P.vcFile_prm);

    set(0, 'UserData', S0);
    save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
    delete_(strrep(P.vcFile_prm, '.prm', '_log.mat')); %delete log file when detecting
end %func
