function [pc1, pc2, pc3] = pcFeaturesBySpike(obj, spikes, sites)
    %PCFEATURESBYSPIKE Get template PC features for selected spikes
    [pc2, pc3] = deal([]);

    nSites = numel(sites);
    nSpikes = numel(spikes);

    pc1 = zeros(nSites, nSpikes);
    if nargout > 1
        pc2 = zeros(nSites, nSpikes);
    elseif nargout > 2
        pc3 = zeros(nSites, nSpikes);
    end
    iFeatures = obj.pcFeatures(:, :, spikes);
    uniqueTemplates = unique(obj.spikeTemplates(spikes));

    % for each template, find the intersection of sites of interest and sites for which we have features
    iSites = ismember(obj.pcFeatureInd(:, uniqueTemplates), sites);
    for iiTemplate = 1:numel(uniqueTemplates)
        iTemplate = uniqueTemplates(iiTemplate);
        spikeIndices = obj.spikeTemplates(spikes) == iTemplate;

        iTempSites = find(iSites(:, iiTemplate));
        tempSites = obj.pcFeatureInd(iTempSites, iTemplate);

        siteIndices = ismember(sites, tempSites);        
        [~, sortIndices] = sort(tempSites);

        pc1(siteIndices, spikeIndices) = iFeatures(1, iTempSites(sortIndices), spikeIndices);
        pc1(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        
        if nargout > 1
            pc2(siteIndices, spikeIndices) = iFeatures(2, iTempSites(sortIndices), spikeIndices);
            pc2(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        end

        if nargout > 2
            pc3(siteIndices, spikeIndices) = iFeatures(3, iTempSites(sortIndices), spikeIndices);
            pc3(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        end
    end
end

