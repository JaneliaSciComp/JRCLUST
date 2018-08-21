%--------------------------------------------------------------------------
function S_clu = S_clu_select_(S_clu, viKeep_clu)
    % automatically trim clusters
    % 7/20/17 JJJ: auto selecting vectors and matrics
    % excl vnSpk_clu, clusterSites, clusterXPositions, clusterYPositions

    % Quality
    fieldNames = fieldnames(S_clu);
    if isempty(fieldNames)
        return;
    end

    viMatch_v = cellfun(@(vi) ~isempty(vi), cellfun(@(cs) regexp(cs, '^v\w*_clu$'), fieldNames, 'UniformOutput', false));
    S_clu = subsetStructFields(S_clu, fieldNames(viMatch_v), viKeep_clu);

    viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), fieldNames, 'UniformOutput', false));
    S_clu = subsetStructFields(S_clu, fieldNames(viMatch_t), viKeep_clu, 3);

    viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), fieldNames, 'UniformOutput', false));
    S_clu = subsetStructFields(S_clu, fieldNames(viMatch_c), viKeep_clu);

    % remap mrWavCor
    if isfield(S_clu, 'mrWavCor')
        S_clu.mrWavCor = S_clu_wavcor_remap_(S_clu, viKeep_clu);
    end
end %func
