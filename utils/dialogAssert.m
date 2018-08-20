%--------------------------------------------------------------------------
% 10/18/17 JJJ: Created
function condition = dialogAssert(condition, msg)
    if ~condition
        errordlg(msg);
        disperr_(msg, '');
    end
end
