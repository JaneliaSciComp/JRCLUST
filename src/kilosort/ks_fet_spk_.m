function [pc1, pc2, pc3] = ks_fet_spk_(spikes, sites, S0)
    %KS_FET_SPK_ get Kilosort features for spikes

    [pc1, pc2, pc3] = deal([]);

    if nargin < 3 || isempty(S0)
        S0 = get0_();
    end
    
    P = S0.P;

    if ~get_set_(P, 'fImportKsort', 0)
        return;
    end

    rez = S0.rez;
    S_clu = S0.S_clu;

    nSites = numel(sites);
    nSpikes = numel(spikes);

    pc1 = zeros(nSites, nSpikes);
    if nargout > 1
        pc2 = zeros(nSites, nSpikes);
    elseif nargout > 2
        pc3 = zeros(nSites, nSpikes);
    end

    pcFeatureIndices = rez.iNeighPC;
    spikeFeatures = permute(rez.cProjPC(spikes, :, :), [2, 3, 1]); % nPCS x nSites x nSpikes
    spikeTemplates = unique(S_clu.viTemplate_spk(spikes));

    % for each template, find the intersection of sites of interest and sites for which we have features
    iSites = ismember(pcFeatureIndices(:, spikeTemplates), sites);
    for iTemp = 1:numel(spikeTemplates)
        template = spikeTemplates(iTemp);
        spikeIndices = S_clu.viTemplate_spk(spikes) == template;

        iTempSites = find(iSites(:, iTemp));
        tempSites = pcFeatureIndices(iTempSites, template);

        siteIndices = ismember(sites, tempSites);        
        [~, sortIndices] = sort(tempSites);

        pc1(siteIndices, spikeIndices) = spikeFeatures(1, iTempSites(sortIndices), spikeIndices);
        pc1(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        
        if nargout > 1
            pc2(siteIndices, spikeIndices) = spikeFeatures(2, iTempSites(sortIndices), spikeIndices);
            pc2(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        end

        if nargout > 2
            pc3(siteIndices, spikeIndices) = spikeFeatures(3, iTempSites(sortIndices), spikeIndices);
            pc3(~siteIndices, spikeIndices) = nan; % PCs not occurring on this site get NaN
        end
    end
end
