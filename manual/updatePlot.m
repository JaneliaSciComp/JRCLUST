%--------------------------------------------------------------------------
function updatePlot(hPlot, newX, newY, userData, tol)
    % update hPlot with new x and y values

    if nargin < 4
        userData = [];
    end
    if nargin < 5
        tol = 1e-10;
    end

    if isempty(hPlot)
        return;
    end

    % set XData and YData to NaN (invisible)
    if isempty(newY) || isempty(newX)
        clearPlots(hPlot);
        return;
    end

    oldX = get(hPlot, 'XData');
    oldY = get(hPlot, 'YData');

    doUpdate = 1;
    if (numel(oldX) == numel(newX)) && (numel(oldY) == numel(newY))
        % don't update if both sets of values haven't changed
        if (norm(oldX(:) - newX(:)) < tol) && (norm(oldY(:) - newY(:)) < tol)
            doUpdate = 0;
        end
    end

    if doUpdate
        set(hPlot, 'XData', newX, 'YData', newY);
    end

    if ~isempty(userData)
        set(hPlot, 'UserData', userData);
    end
end %func
