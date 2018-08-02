%--------------------------------------------------------------------------
function S_clu = fet2clu_(S0, P)
    % can process different shanks separately
    fprintf('Clustering\n');
    S_clu = cluster_spacetime_(S0, P);
    S_clu = postCluster_(S_clu, P);

    if get_set_(P, 'fRepeat_clu', 0), S_clu = S_clu_reclust_(S_clu, S0, P); end

    S_clu = post_merge_(S_clu, P, 0);

    S_clu.spikeClustersAuto = S_clu.spikeClusters;
    fprintf('\tClustering took %0.1f s\n', S_clu.t_runtime);
end %func
