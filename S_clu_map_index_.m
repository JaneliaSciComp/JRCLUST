%--------------------------------------------------------------------------
function S_clu = S_clu_map_index_(S_clu, viMap_clu)
    % update viClu
    vlPos = S_clu.viClu > 0;
    viMap_clu = int32(viMap_clu);
    S_clu.viClu(vlPos) = viMap_clu(S_clu.viClu(vlPos)); %translate cluster number
    % S_clu = S_clu_refresh_(S_clu, 0); % computational efficiency
    % S_clu = S_clu_count_(S_clu);
    S_clu.cviSpk_clu = arrayfun(@(iClu)find(S_clu.viClu==iClu), 1:S_clu.nClu, 'UniformOutput', 0);
    S_clu.vnSpk_clu = cellfun(@numel, S_clu.cviSpk_clu);
    viSite_spk = get0_('viSite_spk');
    S_clu.viSite_clu = double(arrayfun(@(iClu)mode(viSite_spk(S_clu.cviSpk_clu{iClu})), 1:S_clu.nClu));
end %func
