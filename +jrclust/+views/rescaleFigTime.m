function rescaleFigTime(hFigTime, timeScale)
    %RESCALEFIGTIME
    YLim = hFigTime.axApply('default', @get, 'YLim');
    hFigTime.axApply('default', @set, 'YLim', [0, YLim(2)*timeScale]);
    imrectSetPosition(hFigTime, 'hRect', [], [0, YLim(2)*timeScale]);
end
