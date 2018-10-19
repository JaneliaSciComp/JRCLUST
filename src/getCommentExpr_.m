%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function vcComment = getCommentExpr_(vcExpr)
    % Return the comment part of the Matlab code

    iStart = strfind(vcExpr, '%');
    if isempty(iStart), vcComment = ''; return; end
    vcComment = vcExpr(iStart(1):end);
end %func
