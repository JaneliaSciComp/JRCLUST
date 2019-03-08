function fVis = toggleVisible(obj, plotKey)
    %TOGGLEVISIBLE Toggle visibility of plot by key
    if isempty(plotKey)
        return;
    end

    if iscell(plotKey)
        cellfun(@(pKey) obj.toggleVisible(pKey), plotKey);
        return;
    end

    if ~obj.hasPlot(plotKey)
        return;
    end

    hPlot = obj.hPlots(plotKey);
    fVis = jrclust.utils.ifEq(strcmp(get(hPlot, 'Visible'), 'on'), 0, 1);

    % keep this plot hidden (or visible) even if dragging
    if fVis
        obj.permaHidden = setdiff(obj.permaHidden, plotKey);
    else
        obj.permaHidden = union(obj.permaHidden, plotKey);
    end

    obj.setVisible(plotKey, fVis);
end