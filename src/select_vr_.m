%--------------------------------------------------------------------------
function varargout = select_vr_(varargin)
    % [var1, var2, ...] = select_vr(var1, var2, ..., index)

    % sort ascend
    viKeep = varargin{end};
    if islogical(viKeep), viKeep = find(viKeep); end
    for i=1:(nargin-1)
        if isvector(varargin{i})
            varargout{i} = varargin{i}(viKeep);
        else
            varargout{i} = varargin{i}(viKeep, :);
        end
    end
end %func
