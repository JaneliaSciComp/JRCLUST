%--------------------------------------------------------------------------
function [S_clu, nClu_merged] = S_clu_wavcor_merge_(S_clu, P)

    mrWavCor = S_clu.mrWavCor;
    nClu = size(mrWavCor, 2);

    % Identify clusters to remove, update and same (no change), disjoint sets
    % mrWavCor(tril(true(nClu)) | mrWavCor==0) = nan; %ignore bottom half
    mrWavCor(tril(true(nClu))) = 0; %ignore bottom half
    [vrCor_max, viMap_clu] = max(mrWavCor);
    vlKeep_clu = vrCor_max < get_set_(P, 'maxWavCor', 1); % | isnan(vrCor_max);
    if all(vlKeep_clu), nClu_merged=0; return ;end
    min_cor = min(vrCor_max(~vlKeep_clu));
    viClu_same = find(vlKeep_clu);
    viMap_clu(vlKeep_clu) = viClu_same;
    viClu_same = setdiff(viClu_same, viMap_clu(~vlKeep_clu));
    viClu_remove = setdiff(1:nClu, viMap_clu);
    viClu_update = setdiff(setdiff(1:nClu, viClu_same), viClu_remove);
    % viClu_update = setdiff(1:nClu, viClu_same);

    % update cluster number
    try S_clu.icl(viClu_remove) = []; catch, end
    S_clu = S_clu_map_index_(S_clu, viMap_clu); %index mapped
    P.fVerbose = 0;
    S_clu = S_clu_refrac_(S_clu, P); % remove refrac spikes

    % update cluster waveforms and distance
    S_clu = S_clu_wav_(S_clu, viClu_update); %update cluster waveforms
    S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update);
    S_clu = S_clu_remove_empty_(S_clu);

    nClu_merged = nClu - S_clu.nClu;
    fprintf('\tnClu: %d->%d (%d merged, min-cor: %0.4f)\n', nClu, S_clu.nClu, nClu_merged, min_cor);
end %func
