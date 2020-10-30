function success = deleteFile(filename, msg)
%DELETEFILE Delete a file with confirmation.
success = 1;

if nargin < 2
    dlgAns = 'Yes';
else
    dlgAns = questdlg(msg, 'Confirm File Deletion', 'Yes');
end

if strcmp(dlgAns, 'Yes')
    try
        delete(filename);
    catch
        success = 0;
    end
end
end %fun
