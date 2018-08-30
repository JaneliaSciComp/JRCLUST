%--------------------------------------------------------------------------
function S_clu = S_clu_select_(S_clu, viKeep_clu)
    % automatically trim clusters
    % 7/20/17 JJJ: auto selecting vectors and matrices
    % excl nSpikesPerCluster, clusterSites, clusterXPositions, clusterYPositions

    % Quality
    % fieldNames = fieldnames(S_clu);
    % if isempty(fieldNames)
    %     return;
    % end

    % single-dimension reordering
    fieldNames1Dim = clusterFieldsByDim(1);

    S_clu = subsetStructElements(S_clu, fieldNames1Dim, viKeep_clu);

    fieldNames3Dim = clusterFieldsByDim(3);
    S_clu = subsetStructElements(S_clu, fieldNames3Dim, viKeep_clu, 3);

    % remap mrWavCor
    if isfield(S_clu, 'mrWavCor')
        S_clu.mrWavCor = S_clu_wavcor_remap_(S_clu, viKeep_clu);
    end
end % function
