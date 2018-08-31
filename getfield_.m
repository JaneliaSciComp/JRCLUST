%--------------------------------------------------------------------------
function varargout = getfield_(S, varargin)
    for iField = 1:numel(varargin)
        if isfield(S, varargin{iField})
            varargout{iField} = getfield(S, varargin{iField});
        else
            varargout{iField} = [];
        end
    end
end % function
