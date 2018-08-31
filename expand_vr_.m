%--------------------------------------------------------------------------
% 8/16/17 JJJ: created and tested
function vr1 = expand_vr_(vr, nwin, dimm1);
    if nargin<3, dimm1 = [numel(vr) * nwin, 1]; end
    if islogical(vr)
        vr1 = false(dimm1);
    else
        vr1 = zeros(dimm1, 'like', vr);
    end
    vr = repmat(vr(:)', [nwin, 1]);
    vr = vr(:);
    [n,n1] = deal(numel(vr), numel(vr1));
    if n1 > n
        vr1(1:n) = vr;
        vr1(n+1:end) = vr1(n);
    elseif numel(vr1) < n
        vr1 = vr(1:n1);
    else
        vr1 = vr;
    end
end % function
