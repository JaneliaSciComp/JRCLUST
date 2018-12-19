%function [mrFet1, viSpk1] = getDispFeaturesCluster(iClu1, iSite, S0)
function [sampledFeatures, sampledSpikes] = getDispFeaturesCluster(hClust, iCluster, iSite)
    %GETDISPFEATURESCLUSTER Get 
    MAX_SAMPLE = 10000; % max points to display

    hCfg = hClust.hCfg;
    spikeSites = hClust.spikeSites;

    if isempty(iCluster) % select spikes based on sites
        nSites = 1 + round(hCfg.nSiteDir);
        neighbors = hCfg.siteNeighbors(1:nSites, iSite);

        if isempty(hClust.spikesBySite)
            sampledSpikes = find(ismember(spikeSites, neighbors));
        else
            sampledSpikes = jrclust.utils.neCell2mat(hClust.spikesBySite(neighbors)');            
        end

        sampledSpikes = jrclust.utils.subsample(sampledSpikes, MAX_SAMPLE);
    else % get all sites from the cluster
        sampledSpikes = hClust.spikesByCluster{iCluster};
    end

    if strcmp(hCfg.dispFeature, 'vpp')
        sampledWaveforms = squeeze(jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust, sampledSpikes, iSite), hCfg));
        sampledFeatures = max(sampledWaveforms) - min(sampledWaveforms);
    elseif strcmp(hCfg.dispFeature, 'cov')
        sampledFeatures = getSpikeCov(hClust, sampledSpikes, iSite);
    elseif strcmp(hCfg.dispFeature, 'pca')
        sampledFeatures = pca_pc_spk_(sampledSpikes, iSite);
    elseif strcmp(hCfg.dispFeature, 'ppca')
        [mrPv1, mrPv2] = pca_pv_clu_(iSite, S0.iCluCopy);
        sampledFeatures = pca_pc_spk_(sampledSpikes, iSite, mrPv1, mrPv2);
    elseif strcmp(hCfg.dispFeature, 'kilosort')
        sampledFeatures = ks_fet_spk_(sampledSpikes, iSite, S0);
    else
        error('not implemented yet');
    end

    sampledFeatures = squeeze(abs(sampledFeatures));
end
