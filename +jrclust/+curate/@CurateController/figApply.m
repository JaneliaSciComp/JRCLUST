function res = figApply(obj, hFun, varargin)
    %FIGAPPLY Apply a function to all figures
    if nargout == 0 % "Too many output arguments"
        cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs), varargin{:});
    else
        res = cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs), varargin{:});
    end
end