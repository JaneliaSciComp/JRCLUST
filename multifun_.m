%--------------------------------------------------------------------------
function varargout = multifun_(hFun, varargin)
    % apply same function to the input, unary function only

    if nargout ~= numel(varargin), error('n arg mismatch'); end
    for i=1:nargout
        try
            varargout{i} = hFun(varargin{i});
        catch
            varargout{i} = varargin{i};
        end
    end
end %func
