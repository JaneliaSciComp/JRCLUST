%--------------------------------------------------------------------------
% like gcf but doesn't create a new figure
function hFig = gcf_()
    hFig = get(groot(),'CurrentFigure');
end %func
