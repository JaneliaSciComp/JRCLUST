%--------------------------------------------------------------------------
function [spikeTimes11, spikeSites11] = filter_spikes_(spikeTimes0, spikeSites0, tlim)
    % Filter spikes that is within tlim specified

    [spikeTimes11, spikeSites11] = deal([]);
    if isempty(spikeTimes0), return; end
    viKeep11 = find(spikeTimes0 >= tlim(1) & spikeTimes0 <= tlim(end));
    spikeTimes11 = spikeTimes0(viKeep11)  + (1 - tlim(1)); % shift spike timing
    if ~isempty(spikeSites0)
        spikeSites11 = spikeSites0(viKeep11);
    else
        spikeSites11 = [];
    end
end %func
