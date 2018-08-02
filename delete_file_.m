%--------------------------------------------------------------------------
function flag = delete_file_(vcFile)
    flag = 0;
    if ~fileExists(vcFile), return; end
    try
        delete(vcFile);
    catch
        disperr_();
        return;
    end
    % vcCmd = ifeq_(ispc(), 'del', 'rm');
    % eval(sprintf('system(''%s "%s"'');', vcCmd, vcFile));
    flag = 1;
end %func
