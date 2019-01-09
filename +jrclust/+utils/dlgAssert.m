function pred = dlgAssert(pred, failMsg)
    %DLGASSERT Assert with a dialog failure
    if ~pred
        errordlg(failMsg);
        error(failMsg, '');
    end
end
