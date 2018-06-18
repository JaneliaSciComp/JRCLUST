%--------------------------------------------------------------------------
function S_clu = S_clu_combine_(S_clu, S_clu_redo, vlRedo_clu, vlRedo_spk)
    viSpk_cluA = find(~vlRedo_spk);
    [~, viCluA] = ismember(S_clu.viClu(viSpk_cluA), find(~vlRedo_clu));
    S_clu.viClu(viSpk_cluA) = viCluA;
    S_clu.viClu(vlRedo_spk) = S_clu_redo.viClu + max(viCluA);
end %func
