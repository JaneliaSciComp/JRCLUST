function [spikeSites2, spikeSites3] = findSecondaryPeaks(obj, spikeWindows, spikeSites)
    %FINDSECONDARYPEAKS Find secondary peaks
    %   spikeWindows: nSamples x nSites x nSpikes
    evtNeighbors = 2:obj.hCfg.nSitesEvt; % sites within merge radius, excluding ref sites
    siteNeighbors2 = obj.hCfg.siteNeighbors(evtNeighbors, spikeSites);

    spikeWindows2 = spikeWindows(:, evtNeighbors, :);
    spikeVpp = squeeze(max(spikeWindows2) - min(spikeWindows2)); % use Vpp to determine second peak site

    % happens with obj.hCfg.nSitesEvt == 2
    if iscolumn(spikeVpp)
        spikeVpp = spikeVpp';
        peakiestSites = ones(size(spikeVpp), 'like', spikeVpp);
    else
        [~, peakiestSites] = sort(spikeVpp, 'descend');
    end

    spikeSites2 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(1, :), []));

    if size(spikeVpp, 1) > 1
        spikeSites3 = int32(jrclust.utils.rowColSelect(siteNeighbors2, peakiestSites(2, :), []));
    else
        spikeSites3 = [];
    end
end