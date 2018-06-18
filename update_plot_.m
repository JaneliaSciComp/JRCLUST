%--------------------------------------------------------------------------
function update_plot_(hPlot, vrX, vrY, S_plot)
    % update the plot with new x and y

    if nargin<4, S_plot = []; end
    if isempty(hPlot), return; end
    % selective plot to speed up plotting speed
    if isempty(vrY) || isempty(vrX)
        hide_plot_(hPlot);
        %     set(hPlot, 'XData', nan, 'YData', nan);  %visible off
        return;
    end

    % only update if both x and y are changed
    vrX1 = get(hPlot, 'XData');
    vrY1 = get(hPlot, 'YData');
    fUpdate = 1;
    if (numel(vrX1) == numel(vrX)) && (numel(vrY1) == numel(vrY))
        if (std(vrX1(:) - vrX(:)) == 0) && (std(vrY1(:) - vrY(:)) == 0)
            fUpdate = 0;
        end
    end
    if fUpdate, set(hPlot, 'xdata', vrX, 'ydata', vrY); end
    if ~isempty(S_plot), set(hPlot, 'UserData', S_plot); end
end %func
