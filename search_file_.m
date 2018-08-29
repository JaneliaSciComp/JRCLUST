%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile = search_file_(vcFile, csDir)
    % Search file in the provided directory location if it doesn't exist
    if fileExists(vcFile), return ;end

    if ischar(csDir), csDir = {csDir}; end
    for iDir = 1:numel(csDir)
        vcFile_ = replaceDir(vcFile, csDir{iDir});
        if fileExists(vcFile_)
            vcFile = vcFile_;
            return;
        end
    end
    vcFile = []; % file not found
end % function
