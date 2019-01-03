%--------------------------------------------------------------------------
function flag = delete_file_(vcFile)
    flag = 0;
    if ~exist_file_(vcFile), return; end
    try
        delete(vcFile);
    catch
        disperr_();
        return;
    end
    % vcCmd = jrclust.utils.ifEq(ispc(), 'del', 'rm');
    % eval(sprintf('system(''%s "%s"'');', vcCmd, vcFile));
    flag = 1;
end %func
