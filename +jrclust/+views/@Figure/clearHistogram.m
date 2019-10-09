function clearHistogram(obj, plotKey)
    %CLEARPLOT Set XData and YData of a plot to nan
    obj.updateHistogram(plotKey, nan); % updatePlot checks for existence of plotKey
end