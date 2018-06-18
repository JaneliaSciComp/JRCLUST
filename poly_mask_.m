%--------------------------------------------------------------------------
function vlIn = poly_mask_(hPoly, vxPlot, vyPlot)
    mrPolyPos = getPosition(hPoly);
    vxPoly = mrPolyPos([1:end,1],1);
    vyPoly = mrPolyPos([1:end,1],2);
    vlIn = inpolygon(vxPlot, vyPlot, vxPoly, vyPoly);
end %func
