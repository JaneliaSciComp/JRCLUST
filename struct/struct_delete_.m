%--------------------------------------------------------------------------
function S = struct_delete_(S, varargin)
    % delete and set to empty

    for i=1:numel(varargin)
        try
            delete(S.(varargin{i}));
            S.(varargin{i}) = [];
        catch
            ;
        end
    end
end %func
