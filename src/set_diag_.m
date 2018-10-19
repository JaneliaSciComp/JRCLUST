%--------------------------------------------------------------------------
% 10/12/17 JJJ: Works for non-square matrix and constant. Tested
function mr = set_diag_(mr, vr)
    n = min(min(size(mr)), numel(vr));
    % m = numel(vr);
    % mr(sub2ind([n,n], 1:m, 1:m)) = vr;
    mr(sub2ind(size(mr), 1:n, 1:n)) = vr(1:n);
end %func
