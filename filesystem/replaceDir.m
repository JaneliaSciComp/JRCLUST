%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_new = replaceDir(filename, newDirname)
    % Substitute dir
    [newDirname, ~, ~] = fileparts(newDirname);
    [vcDir, filename, vcExt] = fileparts(filename);
    vcFile_new = fullfile(newDirname, [filename, vcExt]);
end % function
