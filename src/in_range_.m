%--------------------------------------------------------------------------
function vl = in_range_(vi, cvi)
    vl = false(size(vi));
    if ~iscell(cvi), cvi = {cvi}; end
    for i=1:numel(cvi)
        lim1 = cvi{i};
        vl = vl | (vi >= lim1(1) & vi <= lim1(2));
    end
end %func
