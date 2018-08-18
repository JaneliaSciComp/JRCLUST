%--------------------------------------------------------------------------
function varargout = get_(varargin)
    % retrieve a field. if not exist then return empty
    % [val1, val2] = get_(S, field1, field2, ...)

    if nargin == 0 || isempty(varargin{1})
        varargout{1} = [];
    else
        S = varargin{1};

        for i=2:nargin
            fieldName = varargin{i};
            try
                varargout{i - 1} = S.(fieldName);
            catch
                varargout{i - 1} = [];
            end
        end
    end

end %func
