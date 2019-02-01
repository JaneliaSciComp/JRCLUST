function varargout = tryGather(varargin)
    %TRYGATHER Try to gather gpuArrays
    for i = 1:nargin
        if isa(varargin{i}, 'gpuArray')
            varargout{i} = gather(varargin{i});
            varargin{i} = [];
        else
            varargout{i} = varargin{i};
        end
    end
end
