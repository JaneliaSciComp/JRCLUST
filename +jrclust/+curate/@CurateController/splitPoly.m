function success = splitPoly(obj, hFig, shift)
    %SPLITPOLY Draw a polygon or rectangle around points in a feature
    %projection and split off a new cluster.
    success = 0;
    
    if numel(obj.selected) > 1
        return;
    end

    iCluster = obj.selected;

    if shift
        hFig.addPlot('hPoly', @imrect);

        try
            polyPos = hFig.plotApply('hPoly', @getPosition);
        catch ME
            hFig.rmPlot('hPoly');
            return;
        end

        if isempty(polyPos)
            return;
        end

        xpos = [repmat(polyPos(1), 2, 1); repmat(polyPos(1) + polyPos(3), 2, 1)];
        ypos = [polyPos(2); repmat(polyPos(2) + polyPos(4), 2, 1); polyPos(2)];
        polyPos = [xpos ypos];
    else
        hFig.addPlot('hPoly', @impoly);
        try
            polyPos = hFig.plotApply('hPoly', @getPosition);
        catch ME
            hFig.rmPlot('hPoly');
            return;
        end
    end

    if isempty(polyPos)
        hFig.rmPlot('hPoly');
        return;
    end

    XData = hFig.plotApply('foreground', @get, 'XData');
    YData = hFig.plotApply('foreground', @get, 'YData');
    splitOff = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));

    hFig.rmPlot('hPoly');
    if ~any(splitOff) || all(splitOff)
        return;
    end

    % show the consequences of the proposed split
    hFig.addPlot('hSplit', @line, XData(splitOff), YData(splitOff), ...
        'Color', obj.hCfg.colorMap(3, :), ...
        'Marker', '.', 'LineStyle', 'none');

    dlgAns = questdlg('Split?', 'Confirmation', 'No');

    hFig.rmPlot('hSplit');

    if strcmp(dlgAns, 'Yes')
        if isfield(hFig.figData, 'foreground') % FigProj
            % first rescale features in this polygon if necessary
            s = hFig.figData.initialScale/hFig.figData.boundScale;
            fgXData = hFig.figData.foreground.XData;
            xFloor = floor(fgXData);
            if strcmp(obj.hCfg.dispFeature, 'vpp')
                fgXData = (fgXData - xFloor)*s;
            else
                % remap to [-1, 1] first, then scale, then map back to [0, 1]
                fgXData = jrclust.utils.linmap(fgXData - xFloor, [0, 1], [-1, 1])*s;
                fgXData = jrclust.utils.linmap(fgXData, [-1, 1], [0, 1]);
            end
            fgXData((fgXData <= 0 | fgXData >= 1)) = nan;
            fgXData = fgXData + xFloor;

            fgYData = hFig.figData.foreground.YData;
            yFloor = floor(fgYData);
            if strcmp(obj.hCfg.dispFeature, 'vpp')
                fgYData = (fgYData - yFloor)*s;
            else
                % remap to [-1, 1] first, then scale, then map back to [0, 1]
                fgYData = jrclust.utils.linmap(fgYData - yFloor, [0, 1], [-1, 1])*s;
                fgYData = jrclust.utils.linmap(fgYData, [-1, 1], [0, 1]);
            end
            fgYData((fgYData <= 0 | fgYData >= 1)) = nan;
            fgYData = fgYData + yFloor;
            % get values of ALL foreground features in this polygon
            ii = any(inpolygon(fgXData, fgYData, polyPos(:, 1), polyPos(:, 2)), 2);
            unitPart = {find(ii)};
        else
            unitPart = {find(splitOff)};
        end
        obj.splitCluster(iCluster, unitPart);
    end
    
    success = 1;
end