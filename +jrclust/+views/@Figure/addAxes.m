function hAx = addAxes(obj, axKey, varargin)
    %ADDAXES Create and store new axes for this figure
    if ~obj.isReady
        obj.hFig = figure();
    end

    if obj.hasAxes(axKey)
        return;
    end

    hAx = axes(obj.hFig, varargin{:});
    obj.hAxes(axKey) = hAx;
    obj.axApply(axKey, @hold, 'on');
end