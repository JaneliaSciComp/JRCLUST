%--------------------------------------------------------------------------
% 10/18/17 JJJ: Created
function condition = assert_(condition, msg)
    if ~condition
        errordlg(msg);
        disperr_(msg, '');
    end
end
