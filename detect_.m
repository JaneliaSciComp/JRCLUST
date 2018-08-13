%--------------------------------------------------------------------------
function detect_(P)
    global spikeTraces spikeWaveforms spikeFeatures;
    runtime_detect = tic;
    % Clear memory (S0 is cleared)
    set(0, 'UserData', []);
    [spikeTraces, spikeWaveforms, spikeFeatures] = deal([]);

    S0 = file2spk_(P);

    if get_set_(P, 'fRamCache', 1)
        spikeTraces = load_bin_(strrep(P.paramFile, '.prm', '_traces.bin'), 'int16', S0.traceDims);
        spikeWaveforms = load_bin_(strrep(P.paramFile, '.prm', '_waveforms.bin'), 'int16', S0.waveformDims);
    end
    spikeFeatures = load_bin_(strrep(P.paramFile, '.prm', '_features.bin'), 'single', S0.featureDims);
    S0.mrPos_spk = spk_pos_(S0, spikeFeatures);

    % measure time
    S0.runtime_detect = toc(runtime_detect);
    fprintf('Detection took %0.1fs for %s\n', S0.runtime_detect, P.paramFile);

    set(0, 'UserData', S0);
    save0_(strrep(P.paramFile, '.prm', '_jrc.mat'));
    delete_(strrep(P.paramFile, '.prm', '_log.mat')); %delete log file when detecting
end %func
