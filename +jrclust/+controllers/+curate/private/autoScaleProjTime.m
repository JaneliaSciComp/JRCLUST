function autoScaleProjTime(hClust, hFigProj, hFigTime, selected)
    %AUTOSCALEPROJTIME Automatically scale features in FigProj and FigTime
    autoScaleFigProj(hFigProj, hClust, selected);
    autoScaleFigTime(hFigTime, hClust, selected);
end
