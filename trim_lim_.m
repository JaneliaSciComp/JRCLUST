%--------------------------------------------------------------------------
% find intersection of two limit ranges
function xlim1 = trim_lim_(xlim1, xlim0)
    dx = diff(xlim1);

    if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
    if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
    xlim1(1) = max(xlim1(1), xlim0(1));
    xlim1(2) = min(xlim1(2), xlim0(2));
end % function
