%--------------------------------------------------------------------------
function set_axis_(hFig, xlim1, ylim1, xlim0, ylim0)
    % set the window within the box limit
    if nargin <= 3
        % square case
        xlim0 = ylim1;
        ylim1 = xlim1;
        ylim0 = xlim0;
    end

    hFig_prev = gcf;
    figure(hFig);
    dx = diff(xlim1);
    dy = diff(ylim1);

    if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
    if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
    if ylim1(1)<ylim0(1), ylim1 = ylim0(1) + [0, dy]; end
    if ylim1(2)>ylim0(2), ylim1 = ylim0(2) + [-dy, 0]; end

    xlim1(1) = max(xlim1(1), xlim0(1));
    ylim1(1) = max(ylim1(1), ylim0(1));
    xlim1(2) = min(xlim1(2), xlim0(2));
    ylim1(2) = min(ylim1(2), ylim0(2));

    axis_([xlim1, ylim1]);

    figure(hFig_prev);
end % function
