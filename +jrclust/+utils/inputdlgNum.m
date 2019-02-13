function dlgAns = inputdlgNum(prompt, dlgTitle, defAns)
    %INPUTDLGNUM Numeric input dialog
    dlgAns = inputdlg(prompt, dlgTitle, 1, {num2str(defAns)});
    try
        dlgAns = str2double(dlgAns{1});
    catch
        dlgAns = nan;
    end
end
