function fVis = setVisible(obj, plotKey, fVis)
    %SETVISIBLE Set visibility of plot by key
    if isempty(plotKey) || nargin < 3
        return;
    end

    if iscell(plotKey)
        cellfun(@(pKey) obj.setVisible(pKey, fVis), plotKey);
        return;
    end

    if ~obj.hasPlot(plotKey)
        return;
    end

    hPlot = obj.hPlots(plotKey);
    set(hPlot, 'Visible', jrclust.utils.ifEq(fVis, 'on', 'off'));
end