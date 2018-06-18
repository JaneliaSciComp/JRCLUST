%--------------------------------------------------------------------------
function S = select_polygon_(hPlot)
    S = get(hPlot, 'UserData');
    % try delete(S.hPlot_split); catch; end
    update_plot2_proj_();
    % try delete(S.hPoly); catch; end

    S.hPoly = impoly_(); %get a polygon drawing from user
    if isempty(S.hPoly), S=[]; return; end;
    mrPolyPos = getPosition(S.hPoly);

    vrXp = get(hPlot, 'XData');
    vrYp = get(hPlot, 'YData');
    vlKeep1 = inpolygon(vrXp, vrYp, mrPolyPos(:,1), mrPolyPos(:,2));
    % viKeep1 = find(inpoly([vrXp(:), vrYp(:)], mrPolyPos));

    [viEvtPlot,~,~] = ind2sub(S.tr_dim, S.viPlot);
    viEvtKeep1 = unique(viEvtPlot(vlKeep1));
    viEvtKeep = find(ismember(viEvtPlot, viEvtKeep1));

    update_plot2_proj_(vrXp(viEvtKeep), vrYp(viEvtKeep));
    % S.hPlot_split = line(vrXp(viEvtKeep), vrYp(viEvtKeep), 'Color', [1 0 0], 'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none');
    % nSites = S.tr_dim(2);
    % axis([0 nSites 0 nSites]); %reset view

    set(hPlot, 'UserData', S);
    hold off;
end %fund
