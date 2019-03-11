function updateFigTime(obj, doAutoscale)
    %UPDATEFIGTIME
    if ~obj.hasFig('FigTime')
        return;
    end

    hFigTime = obj.hFigs('FigTime');
    jrclust.views.plotFigTime(hFigTime, obj.hClust, obj.hCfg, obj.selected, obj.maxAmp, obj.currentSite);
    hFigTime.setMouseable(); % no special mouse function

    if doAutoscale
        jrclust.views.autoScaleFigTime(hFigTime, obj.hClust, obj.selected);
    end
end