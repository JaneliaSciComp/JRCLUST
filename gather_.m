%--------------------------------------------------------------------------
function varargout = gather_(varargin)
    for i=1:nargin
        varargout{i} = varargin{i};
        if isa(varargin{i}, 'gpuArray')
            try
                varargout{i} = gather(varargin{i});
            catch
                ;
            end
        end
    end
end %func
