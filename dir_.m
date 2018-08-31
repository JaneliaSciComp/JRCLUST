%--------------------------------------------------------------------------
function [csFiles_full, csFiles] = dir_(vcFilter_dir, csExcl)
    % return name of files full path, exclude files
    if nargin>=2
        if ischar(csExcl), csExcl = {csExcl}; end
        csExcl = union(csExcl, {'.', '..'});
    else
        csExcl = [];
    end
    csFiles = dir(vcFilter_dir);
    csFiles  = {csFiles.('name')};
    csFiles = setdiff(csFiles, csExcl);
    [vcDir, ~, ~] = fileparts(vcFilter_dir);
    if isempty(vcDir), vcDir='.'; end
    csFiles_full = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end % function
