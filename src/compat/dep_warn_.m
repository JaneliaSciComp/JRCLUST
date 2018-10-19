%--------------------------------------------------------------------------
function dep_warn_(cmd, helpStr)

    fprintf(2, '`%s` has been deprecated\n', cmd);
    if nargin > 1 && ~isempty(helpStr)
        fprintf(2, '%s\n', helpStr);
    end

end % func
