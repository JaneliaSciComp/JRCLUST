%--------------------------------------------------------------------------
function [vr, vi] = get_diag_(mr)
    n = min(size(mr));
    vi = sub2ind([n,n], 1:n, 1:n);
    vr = mr(vi);
end %func
