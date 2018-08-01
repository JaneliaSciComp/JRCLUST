%--------------------------------------------------------------------------
function hFig = resizeFigure(hFig, initialPos, refocus)
    if nargin < 3
        refocus = 1;
    end

    taskbarHeight = 40;

    pos0 = get(groot, 'ScreenSize');
    width = pos0(3);
    height = pos0(4) - taskbarHeight;

    finalPos = [0 0 0 0];
    finalPos(1) = max(round(initialPos(1) * width), 1);
    finalPos(2) = max(round(initialPos(2) * height), 1) + taskbarHeight;
    finalPos(3) = min(round(initialPos(3) * width), width);
    finalPos(4) = min(round(initialPos(4) * height), height);

    if isempty(hFig)
        hFig = figure; % create figure
    else
        hFig = figure(hFig);
    end

    drawnow;

    set(hFig, 'OuterPosition', finalPos, 'Color', 'w', 'NumberTitle', 'off');
end
