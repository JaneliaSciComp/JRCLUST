%--------------------------------------------------------------------------
function S_clu = S_clu_update_note_(S_clu, iClu1, iClu_next)
    if isempty(iClu_next), return ;end
    vcNote_clu1 = S_clu.csNote_clu{iClu1};
    if isempty(vcNote_clu1), return; end
    iStart = find(vcNote_clu1 == '=', 1, 'first');
    if isempty(iStart), return; end
    vcPre = vcNote_clu1(1:iStart);
    vcNote_clu1 = vcNote_clu1(iStart+1:end);
    iEnd = find(vcNote_clu1 == ',' | vcNote_clu1 == ';' | vcNote_clu1 == ' ', 1, 'first');
    if ~isempty(iEnd)
        vcNote_clu1 = vcNote_clu1(1:iEnd-1);
        vcPost = vcNote_clu1(iEnd:end);
    else
        vcPost = '';
    end
    S_clu.csNote_clu{iClu1} = sprintf('%s%d%s', vcPre, iClu_next, vcPost);

    if isnan(iClu_next), iClu_next = []; return; end
end %func
