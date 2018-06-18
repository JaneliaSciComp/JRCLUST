%--------------------------------------------------------------------------
% 10/18/17 JJJ: Created
function flag = assert_(flag, vcMsg)
    if ~flag
        errordlg(vcMsg);
        disperr_(vcMsg, '');
    end
end
