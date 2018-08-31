%--------------------------------------------------------------------------
function plot_update_(vhPlot, mrX, mrY)
    nPlots = numel(vhPlot);
    for iPlot=1:nPlots
        vrX1 = mrX(:, iPlot:nPlots:end);
        vrY1 = mrY(:, iPlot:nPlots:end);
        set(vhPlot(iPlot), 'XData', vrX1(:), 'YData', vrY1(:));
    end
end % function
