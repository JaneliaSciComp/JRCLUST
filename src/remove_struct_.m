%--------------------------------------------------------------------------
function S = remove_struct_(S, varargin)
    % remove fields from a struct
    for i=1:numel(varargin)
        if isfield(S, varargin{i})
            S = rmfield(S, varargin{i});
        end
    end
end %func
