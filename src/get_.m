%--------------------------------------------------------------------------
function varargout = get_(varargin)
    % retrieve a field. if not exist then return empty
    % [val1, val2] = get_(S, field1, field2, ...)

    if nargin==0, varargout{1} = []; return; end
    S = varargin{1};
    if isempty(S), varargout{1} = []; return; end

    for i=2:nargin
        vcField = varargin{i};
        try
            varargout{i-1} = S.(vcField);
        catch
            varargout{i-1} = [];
        end
    end
end %func
