%--------------------------------------------------------------------------
function S = struct_(varargin)
    % smart about dealing with cell input
    for iArg = 1:2:numel(varargin)
        try
            S.(varargin{iArg}) = varargin{iArg+1};
        catch
            disperr_('struct_');
        end
    end
end
