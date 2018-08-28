%--------------------------------------------------------------------------
function clearPlots(hPlots)
    % set XData and YData to NaN for all input plots

    for i = 1:numel(hPlots)
        set(hPlots(i), 'XData', nan, 'YData', nan);
    end
end % func
