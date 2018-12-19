%function imrect_set_(hRect, xpos, ypos)
function imrect_set_(hFig, plotKey, xpos, ypos)
    %vrPos = getPosition(hRect);
    vrPos = hFig.imrectFun(plotKey, @getPosition);
    if ~isempty(xpos)
        vrPos(1) = min(xpos);
        vrPos(3) = abs(diff(xpos));
    end
    if ~isempty(ypos)
        vrPos(2) = min(ypos);
        vrPos(4) = abs(diff(ypos));
    end
    %setPosition(hRect, vrPos);
    hFig.imrectFun(plotKey, @setPosition, vrPos);
end
