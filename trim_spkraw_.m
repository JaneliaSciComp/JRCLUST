%--------------------------------------------------------------------------
function tmrWav = trim_spkraw_(tmrWav, P)
    spkLim_raw = get_(P, 'spkLim_raw');
    nSamples_raw = diff(spkLim_raw) + 1;
    spkLim_factor_merge = getOr(P, 'spkLim_factor_merge', 1);
    spkLim_merge = round(P.spkLim * spkLim_factor_merge);
    nSamples_raw_merge = diff(spkLim_merge) + 1;
    if size(tmrWav,1) <= nSamples_raw_merge, return ;end

    lim_merge = [spkLim_merge(1) - spkLim_raw(1) + 1,  nSamples_raw - spkLim_raw(2) + spkLim_merge(2)];
    tmrWav = tmrWav(lim_merge(1):lim_merge(2), :, :);
    tmrWav = meanSubt_(tmrWav);
end %func
