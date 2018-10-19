%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile = search_file_(vcFile, csDir)
    % Search file in the provided directory location if it doesn't exist
    if exist_file_(vcFile)
        return;
    end

    if ischar(csDir)
        csDir = {csDir};
    end

    [dirname, filename, ext] = fileparts(vcFile);
    filename = [filename ext];

    for iDir = 1:numel(csDir)
        dirname = csDir{iDir};
        vcFile_ = fullfile(dirname, filename);

        if exist_file_(vcFile_)
            vcFile = vcFile_;
            return;
        end
    end

    vcFile = []; % file not found
end %func
