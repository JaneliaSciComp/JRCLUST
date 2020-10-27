function success = deleteFile(msg, filename)
%DELETEFILE Delete a file with confirmation.
success = 1;
dlgans = questdlg(msg, 'Confirm File Deletion', 'Yes');

if strcmp(dlgans, 'Yes')
    try
        delete(filename);
    catch
        success = 0;
    end
end
end %fun

