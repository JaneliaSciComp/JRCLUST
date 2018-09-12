%--------------------------------------------------------------------------
function mkdir_(vcDir)
    % make only if it doesn't exist. provide full path for dir
    if ~exist(vcDir, 'dir'), mkdir(vcDir); end
end %func
