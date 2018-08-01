%--------------------------------------------------------------------------
function [spikeTimes11, viSite_spk11] = filter_spikes_(spikeTimes0, viSite_spk0, tlim)
    % Filter spikes that is within tlim specified

    [spikeTimes11, viSite_spk11] = deal([]);
    if isempty(spikeTimes0), return; end
    viKeep11 = find(spikeTimes0 >= tlim(1) & spikeTimes0 <= tlim(end));
    spikeTimes11 = spikeTimes0(viKeep11)  + (1 - tlim(1)); % shift spike timing
    if ~isempty(viSite_spk0)
        viSite_spk11 = viSite_spk0(viKeep11);
    else
        viSite_spk11 = [];
    end
end %func
