%--------------------------------------------------------------------------
function mouse_figure_(hFig, vhAx)
    for iAx = 1:numel(vhAx)
        mouse_figure(hFig, vhAx(iAx));
    end
end %func
