%--------------------------------------------------------------------------
function spkLim_raw = calc_spkLim_raw_(P)
    % Calculate the raw spike waveform range

    spkLim_raw_ms = get_(P, 'spkLim_raw_ms');
    if isempty(spkLim_raw_ms)
        spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
        spkLim_raw = spkLim * get_set_(P, 'spkLim_raw_factor', 2); % backward compatible
    else
        spkLim_raw = round(P.spkLim_raw_ms * P.sRateHz / 1000);
    end
end %func
