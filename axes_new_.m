%--------------------------------------------------------------------------
function hAx = axes_new_(hFig)
    if ischar(hFig), hFig = figureByTag(hFig); end
    figure(hFig); %set focus to figure %might be slow
    clf(hFig);
    hAx = axes();
    hold(hAx, 'on');
end %func
