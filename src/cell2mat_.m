%--------------------------------------------------------------------------
function vr = cell2mat_(cvr)
    % remove empty
    vi = find(cellfun(@(x)~isempty(x), cvr));
    vr = cell2mat(cvr(vi));
end %func
