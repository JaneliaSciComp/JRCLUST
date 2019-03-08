function vals = subplotApply(obj, plotKey, spIndex, hFun, varargin)
    %SUBPLOTAPPLY Apply a plotting function to a subplot
    vals = [];

    if ~obj.hasSubplot(plotKey) || numel(obj.hSubplots(plotKey)) < spIndex
        return;
    end

    hAxes_ = obj.hSubplots(plotKey);
    hAx = hAxes_(spIndex);

    % apply hFun
    if nargout == 1
        vals = hFun(hAx, varargin{:});
    else
        hFun(hAx, varargin{:});
    end

    % save hAxes
    obj.hSubplots(plotKey) = hAxes_;
end