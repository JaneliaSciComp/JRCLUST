%--------------------------------------------------------------------------
function varargout = gather_(varargin)
    for i=1:nargin
        varargout{i} = varargin{i};
        if isa(varargin{i}, 'gpuArray')
            varargout{i} = gather(varargin{i});
        end
    end
end % function
