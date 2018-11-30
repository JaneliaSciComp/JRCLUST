function meanWf = trimRawWaveforms(meanWf, hCfg)
    %TRIMRAWWAVEFORMS Pare raw waveforms down to the size of filtered waveforms
    nSamplesRaw = diff(hCfg.evtWindowRawSamp) + 1;
    spkLimMerge = round(hCfg.evtWindowSamp * hCfg.spkLim_factor_merge);
    nSamplesRawMerge = diff(spkLimMerge) + 1;

    if size(meanWf, 1) <= nSamplesRawMerge
        return;
    end

    limits = [spkLimMerge(1) - hCfg.evtWindowRawSamp(1) + 1, ...
              nSamplesRaw - (hCfg.evtWindowRawSamp(2) - spkLimMerge(2))];

    meanWf = meanWf(limits(1):limits(2), :, :);
    meanWf = jrclust.utils.meanSubtract(meanWf);
end