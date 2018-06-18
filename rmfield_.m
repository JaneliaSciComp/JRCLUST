%--------------------------------------------------------------------------
function S = rmfield_(S, varargin)
    % varargin: list of fields to remove
    for i=1:numel(varargin)
        if isfield(S, varargin{i})
            S = rmfield(S, varargin{i});
        end
    end
end %func
