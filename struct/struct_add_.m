%--------------------------------------------------------------------------
function S = struct_add_(S, varargin)

    for i=1:numel(varargin)
        S.(inputname(i+1)) = varargin{i};
    end
end % function
