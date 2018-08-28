%--------------------------------------------------------------------------
function [featuresMin, featuresMax] = getPCFeatures(spikes, sitesOfInterest, S0)
    rez = S0.rez;
    S_clu = S0.S_clu;
    
    nSites = numel(sitesOfInterest);
    nSpikes = numel(spikes);
    
    featuresMin = zeros(nSites, nSpikes);
    featuresMax = zeros(nSites, nSpikes);

    pcFeatureIndices = rez.iNeighPC;
    spikeFeatures = permute(rez.cProjPC(spikes, :, :), [2 3 1]); % nPCs x nSites x nSpikes
    spikeTemplates = unique(S_clu.spikeTemplates(spikes));

    % for each template, find the intersection of sites of interest and sites for which we have features
    iSites = ismember(pcFeatureIndices(:, spikeTemplates), sitesOfInterest);
    for iTemplate = 1:numel(spikeTemplates)
        template = spikeTemplates(iTemplate);
        spikeIndices = S_clu.spikeTemplates(spikes) == template;

        iTemplateSites = find(iSites(:, iTemplate));
        templateSites = pcFeatureIndices(iTemplateSites, template);
        siteIndices = ismember(sitesOfInterest, templateSites);
        
        [~, sortIndices] = sort(templateSites);
        featuresMin(siteIndices, spikeIndices) = spikeFeatures(S0.kspc(1), iTemplateSites(sortIndices), spikeIndices);
        featuresMin(~siteIndices, spikeIndices) = nan;
        featuresMax(siteIndices, spikeIndices) = spikeFeatures(S0.kspc(2), iTemplateSites(sortIndices), spikeIndices);
        featuresMax(~siteIndices, spikeIndices) = nan;
    end
end % func
