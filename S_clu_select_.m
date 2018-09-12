%--------------------------------------------------------------------------
function S_clu = S_clu_select_(S_clu, viKeep_clu)
    % automatically trim clusters
    % 7/20/17 JJJ: auto selecting vectors and matrics
    % excl vnSpk_clu, viSite_clu, vrPosX_clu, vrPosY_clu

    % Quality
    csNames = fieldnames(S_clu);
    if isempty(csNames), return; end
    viMatch_v = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
    S_clu = struct_select_(S_clu, csNames(viMatch_v), viKeep_clu);

    viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
    S_clu = struct_select_(S_clu, csNames(viMatch_t), viKeep_clu, 3);

    viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
    S_clu = struct_select_(S_clu, csNames(viMatch_c), viKeep_clu);

    % remap mrWavCor
    if isfield(S_clu, 'mrWavCor')
        S_clu.mrWavCor = S_clu_wavcor_remap_(S_clu, viKeep_clu);
    end
end %func
