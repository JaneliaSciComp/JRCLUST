function mouseClickFigWav(obj, xyPos, clickType)
    %MOUSECLICKFIGWAV Handle callbacks for mouse clicks in main view
    if obj.isWorking
        return;
    end

    iCluster = round(xyPos(1)); % floor of x position
    if iCluster < 1 || iCluster > obj.hClust.nClusters
        return;
    end

    if strcmp(clickType, 'normal')
        obj.updateSelect(iCluster);
    elseif strcmp(clickType, 'alt') % right click, select secondary cluster
        obj.updateSelect([obj.selected(1) iCluster]);
    else                            % middle click, ignore
        disp(clickType);
        return;
    end

    hFig = obj.hFigs('FigWav');
    hFig.wait(1);
    hFig.wait(0);
end