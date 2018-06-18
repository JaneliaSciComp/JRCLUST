%--------------------------------------------------------------------------
% 10/27/17: Detailed error report
function flag = S_clu_valid_(S_clu)
    % check the validity of the cluster by making sure the number of clusters is self-consistent.

    flag = 0; %assumne invalid by default
    csNames = fieldnames(S_clu);
    if isempty(csNames), return; end
    viMatch_v = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
    viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
    viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
    viMatch_m = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^m\w*_clu$'), csNames, 'UniformOutput', false));
    [viMatch_v, viMatch_t, viMatch_c, viMatch_m] = multifun_(@find, viMatch_v, viMatch_t, viMatch_c, viMatch_m);
    csNames_m = csNames(viMatch_m);
    csNames_m{end+1} = 'mrWavCor';
    nClu = S_clu.nClu;

    vlError_v = cellfun(@(vc)numel(S_clu.(vc)) ~= nClu, csNames(viMatch_v));
    vlError_t = cellfun(@(vc)size(S_clu.(vc),3) ~= nClu, csNames(viMatch_t));
    vlError_c = cellfun(@(vc)numel(S_clu.(vc)) ~= nClu, csNames(viMatch_c));
    vlError_m = cellfun(@(vc)all(size(S_clu.(vc)) ~= [nClu, nClu]), csNames_m);
    [flag_v, flag_t, flag_c, flag_m] = multifun_(@(x)~any(x), vlError_v, vlError_t, vlError_c, vlError_m);
    flag = flag_v && flag_t && flag_c && flag_m;
    cell2vc__ = @(x)sprintf('%s, ', x{:});
    if ~flag
        fprintf(2, 'S_clu_valid_: flag_v:%d, flag_t:%d, flag_c:%d, flag_m:%d\n', flag_v, flag_t, flag_c, flag_m);
        if ~flag_v, fprintf(2, '\t%s\n', cell2vc__(csNames(viMatch_v(vlError_v)))); end
        if ~flag_t, fprintf(2, '\t%s\n', cell2vc__(csNames(viMatch_t(vlError_t)))); end
        if ~flag_c, fprintf(2, '\t%s\n', cell2vc__(csNames(viMatch_c(vlError_c)))); end
        if ~flag_m, fprintf(2, '\t%s\n', cell2vc__(csNames(viMatch_m(vlError_m)))); end
    end
end %func
