%--------------------------------------------------------------------------
function hide_plot_(vhPlot)
    for i=1:numel(vhPlot), set(vhPlot(i), 'XData', nan, 'YData', nan); end
end %func
