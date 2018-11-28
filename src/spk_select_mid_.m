function subSpikes = spk_select_mid_(spikes, spikeTimes, nTime_clu)
    % viTime_spk = get0_('viTime_spk');
    iSpikeMid = round(numel(spikeTimes)/2); % index of the middlest spike
    nearestToCenter = jrclust.utils.rankorder(abs(spikes - iSpikeMid), 'ascend');
    nSpikesInterval = round(numel(spikes) / nTime_clu);
    subSpikes = spikes(nearestToCenter <= nSpikesInterval);
end
