%--------------------------------------------------------------------------
function varargout = sort_ascend_(varargin)
    % sort all the other fields basedon the first field in ascending order
    % [a', b', c'] = sort_ascend_(a, b, c)

    [varargout{1}, viSrt] = sort(varargin{1}, 'ascend');
    for i=2:nargin
        varargout{i} = varargin{i}(viSrt);
    end
end %func
