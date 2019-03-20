function proj = templateFeaturesBySpike(obj, spikes, sites)
    %TEMPLATEFEATURESBYSPIKE Get template features features for selected spikes on given
    %sites
    nSites = numel(sites);
    nSpikes = numel(spikes);

    proj = zeros(nSites, nSpikes);
    iFeatures = obj.templateFeatures(:, spikes);
    uniqueTemplates = unique(obj.spikeTemplates(spikes));

    % for each template, find the intersection of sites of interest and sites for which we have features
    iSites = ismember(obj.templateFeatureInd(:, uniqueTemplates), sites);
    for iiTemplate = 1:numel(uniqueTemplates)
        iTemplate = uniqueTemplates(iiTemplate);
        spikeIndices = obj.spikeTemplates(spikes) == iTemplate;

        iTempSites = find(iSites(:, iiTemplate));
        tempSites = obj.templateFeatureInd(iTempSites, iTemplate);

        siteIndices = ismember(sites, tempSites);        
        [~, sortIndices] = sort(tempSites);

        proj(siteIndices, spikeIndices) = iFeatures(iTempSites(sortIndices), spikeIndices);
        proj(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
    end
end

