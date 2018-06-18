%--------------------------------------------------------------------------
% 10/24/17 JJJ: remap mrWavCor
function mrWavCor_new = S_clu_wavcor_remap_(S_clu, viKeep_clu)
    if islogical(viKeep_clu), viKeep_clu = find(viKeep_clu); end
    nClu_old = size(S_clu.mrWavCor, 1);
    viOld = find(~isnan(S_clu.mrWavCor));
    [viCol, viRow] = ind2sub([nClu_old,nClu_old], viOld);
    vlKeep = find(ismember(viCol, viKeep_clu) & ismember(viRow, viKeep_clu));
    [viCol, viRow, viOld] = deal(viCol(vlKeep), viRow(vlKeep), viOld(vlKeep));

    nClu_new = numel(viKeep_clu);
    mrWavCor_new = zeros(nClu_new);
    viOld2New = zeros(nClu_old, 1);
    viOld2New(viKeep_clu) = 1:nClu_new;
    viNew = sub2ind([nClu_new, nClu_new], viOld2New(viCol), viOld2New(viRow));
    mrWavCor_new(viNew) = S_clu.mrWavCor(viOld);
end %func
