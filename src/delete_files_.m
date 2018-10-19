%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_files_(csFiles, fVerbose)
    % Delete list of files
    % delete_files_(vcFile)
    % delete_files_(csFiles)
    % delete_files_(csFiles, fVerbose)

    if nargin<2, fVerbose = 1; end
    if ischar(csFiles), csFiles = {csFiles}; end
    for iFile = 1:numel(csFiles)
        fDeleted = delete_file_(csFiles{iFile});
        if fVerbose && fDeleted
            fprintf('\tdeleted %s.\n', csFiles{iFile});
        end
    end
end %func
