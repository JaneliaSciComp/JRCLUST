%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_new = subsDir_(vcFile, vcDir_new)
    % Substitute dir
    [vcDir_new,~,~] = fileparts(vcDir_new);
    [vcDir, vcFile, vcExt] = fileparts(vcFile);
    vcFile_new = fullfile(vcDir_new, [vcFile, vcExt]);
end % func
