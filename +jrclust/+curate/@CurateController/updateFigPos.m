function updateFigPos(obj)
    %UPDATEFIGPOS Plot cluster position on probe
    if isempty(obj.selected) || ~obj.hasFig('FigPos')
        return;
    end

    hFigPos = obj.hFigs('FigPos');
    jrclust.views.plotFigPos(hFigPos, obj.hClust, obj.hCfg, obj.selected, obj.maxAmp);
    hFigPos.hFunKey = @(hO, hE) []; % do-nothing function
    hFigPos.setMouseable();
end