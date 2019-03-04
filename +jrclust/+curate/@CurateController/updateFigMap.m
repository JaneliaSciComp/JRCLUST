function updateFigMap(obj)
    %UPDATEFIGMAP Plot probe map
    if isempty(obj.selected) || ~obj.hasFig('FigMap')
        return;
    end

    hFigMap = obj.hFigs('FigMap');
    jrclust.views.plotFigMap(hFigMap, obj.hClust, obj.hCfg, obj.selected);
    hFigMap.hFunKey = @(hO, hE) []; % do-nothing function
    hFigMap.setMouseable();
end
