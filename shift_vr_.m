%--------------------------------------------------------------------------
function mr = shift_vr_(vr, vn_shift)
    n = numel(vr);
    % mr = zeros(n, numel(vn_shift), 'like', vr);
    mi = bsxfun(@plus, (1:n)', vn_shift(:)');
    mi(mi<1) = 1;
    mi(mi>n) = n;
    mr = vr(mi);
end %func
