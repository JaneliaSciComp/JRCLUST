%--------------------------------------------------------------------------
function S_clu = S_clu_map_index_(S_clu, viMap_clu)
    % update viClu
    vlPos = S_clu.spikeClusters > 0;
    viMap_clu = int32(viMap_clu);
    S_clu.spikeClusters(vlPos) = viMap_clu(S_clu.spikeClusters(vlPos)); %translate cluster number
    % S_clu = S_clu_refresh_(S_clu, 0); % computational efficiency
    % S_clu = S_clu_count_(S_clu);
    S_clu.cviSpk_clu = arrayfun(@(iClu)find(S_clu.spikeClusters==iClu), 1:S_clu.nClusters, 'UniformOutput', 0);
    S_clu.vnSpk_clu = cellfun(@numel, S_clu.cviSpk_clu);
    viSite_spk = get0_('viSite_spk');
    S_clu.clusterSites = double(arrayfun(@(iClu)mode(viSite_spk(S_clu.cviSpk_clu{iClu})), 1:S_clu.nClusters));
end %func
