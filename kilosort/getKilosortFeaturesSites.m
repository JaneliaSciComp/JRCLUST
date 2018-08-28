%--------------------------------------------------------------------------
function [pc1, pc2] = getKilosortFeaturesSites(spikes, sitesOfInterest, S0)
    % return 1st and 2nd (or 1st and 3rd, or 2nd and 3rd) principal components
    % for each spike on sites of interest

    rez = S0.rez;
    S_clu = S0.S_clu;

    nSites = numel(sitesOfInterest);
    nSpikes = numel(spikes);

    pc1 = zeros(nSites, nSpikes);
    pc2 = zeros(nSites, nSpikes);

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
        pc1(siteIndices, spikeIndices) = spikeFeatures(S0.kspc(1), iTemplateSites(sortIndices), spikeIndices);
        pc1(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        pc2(siteIndices, spikeIndices) = spikeFeatures(S0.kspc(2), iTemplateSites(sortIndices), spikeIndices);
        pc2(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
    end
end % func
