%--------------------------------------------------------------------------
function [vr1, vr2] = swap_vr_(vr1, vr2, vi_swap)
    % assert: vr1 and vr2 must have the same dimensions
    [vr1(vi_swap), vr2(vi_swap)] = deal(vr2(vi_swap), vr1(vi_swap));
end
