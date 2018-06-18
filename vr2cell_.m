%--------------------------------------------------------------------------
function cvr = vr2cell_(vr, cvi)
    cvr = cellfun(@(vi)vr(vi), cvi, 'UniformOutput', 0);
end %func
