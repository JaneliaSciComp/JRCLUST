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
    retained = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));

    hFig.rmPlot('hPoly');
    if ~any(retained)
        return;
    end

    % show the consequences of the proposed split
    hFig.addPlot('hSplit', @line, XData(retained), YData(retained), ...
        'Color', obj.hCfg.colorMap(3, :), ...
        'Marker', '.', 'LineStyle', 'none');

    dlgAns = questdlg('Split?', 'Confirmation', 'No');

    hFig.rmPlot('hSplit');

    if strcmp(dlgAns, 'Yes')
        iSpikes = obj.hClust.spikesByCluster{iCluster};
        % convert global (spike-table) indices of split-off spikes to local
        % (within-cluster) indices for unitPart
        if isfield(hFig.figData, 'dispFeatures')
%             unitPart = {find(ismember(iSpikes, hFig.figData.dispFeatures.fgSpikes(retained)))};
            unitPart = {find(ismember(hFig.figData.dispFeatures.fgSpikes(retained), iSpikes))};
        else
            unitPart = {find(retained)};
        end
        obj.splitCluster(iCluster, unitPart);
    end
    
    success = 1;
end