function success = deleteFile(msg, filename)
%DELETEFILE Delete a file with confirmation.
success = 1;
dlgans = questdlg(msg, 'Confirm File Deletion', 'OK', 'Cancel', 'OK');

if strcmp(dlgans, 'OK')
    try
        delete(filename);
    catch
        success = 0;
    end
end
end

