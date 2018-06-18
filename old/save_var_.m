%--------------------------------------------------------------------------
function save_var_(vcFile, varargin)
    % must pass
    struct_save_(struct_(varargin{:}), vcFile);
end
