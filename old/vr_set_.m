%--------------------------------------------------------------------------
function vr = vr_set_(vr, vi, val);
    vi(vi<1 | vi>numel(vr)) = [];
    vr(vi) = val;
end %func
