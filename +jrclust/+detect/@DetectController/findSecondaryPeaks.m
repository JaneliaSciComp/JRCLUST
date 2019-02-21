function [spikeSites2, spikeSites3] = findSecondaryPeaks(obj, spikeWindows, spikeSites)
    %FINDSECONDARYPEAKS Find secondary peaks
    %   spikeWindows: nSamples x nSites x nSpikes
    evtNeighbors = 2:obj.hCfg.nSitesEvt; % sites within merge radius, excluding ref sites
    siteNeighbors2 = obj.hCfg.siteNeighbors(evtNeighbors, spikeSites);

    spikeWindows2 = spikeWindows(:, evtNeighbors, :);
    spikeVpp = squeeze(max(spikeWindows2) - min(spikeWindows2)); % use Vpp to determine second peak site

    [~, peakiestSites] = sort(spikeVpp, 'descend');
    spikeSites2 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(1, :), []));
    spikeSites3 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(2, :), []));
end