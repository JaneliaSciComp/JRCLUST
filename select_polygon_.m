%--------------------------------------------------------------------------
function S = select_polygon_(hPlot)
    S = get(hPlot, 'UserData');

    S.hPoly = impoly_(); %get a polygon drawing from user

    if isempty(S.hPoly)
        S = [];
        return;
    end

    polyPos = getPosition(S.hPoly);

    xvals = get(hPlot, 'XData');
    yvals = get(hPlot, 'YData');
    vlKeep1 = inpolygon(xvals, yvals, polyPos(:,1), polyPos(:,2));

    [viEvtPlot, ~, ~] = ind2sub(S.tr_dim, S.viPlot);
    viEvtKeep1 = unique(viEvtPlot(vlKeep1));
    viEvtKeep = find(ismember(viEvtPlot, viEvtKeep1));

    update_plot2_proj_(xvals(viEvtKeep), yvals(viEvtKeep));

    set(hPlot, 'UserData', S);
    hold off;
end %fund
