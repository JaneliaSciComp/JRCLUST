%--------------------------------------------------------------------------
function hFig = resize_figure_(hFig, posvec0, fRefocus)
    height_taskbar = 40;

    pos0 = get(groot, 'ScreenSize');
    width = pos0(3);
    height = pos0(4) - height_taskbar;
    % width = width;
    % height = height - 132; %width offset
    % width = width - 32;
    posvec = [0 0 0 0];
    posvec(1) = max(round(posvec0(1)*width),1);
    posvec(2) = max(round(posvec0(2)*height),1) + height_taskbar;
    posvec(3) = min(round(posvec0(3)*width), width);
    posvec(4) = min(round(posvec0(4)*height), height);
    % drawnow;
    if isempty(hFig)
        hFig = figure; %create a figure
    else
        hFig = figure(hFig);
    end
    drawnow;
    set(hFig, 'OuterPosition', posvec, 'Color', 'w', 'NumberTitle', 'off');
end
