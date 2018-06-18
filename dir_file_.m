%--------------------------------------------------------------------------
% 8/7/17 JJJ: tested and documented, file creation date
function [csFile_merge, vcDir] = dir_file_(vcFile_dir, fSortByDate)
    % function [csFile_merge, vcDir, csFile_merge1] = dir_file_(vcFile_dir, fSortByDate)
    % search for files and sort by date

    if nargin<2, fSortByDate=1; end

    [vcDir, ~, ~] = fileparts(vcFile_dir);
    if isempty(vcDir), vcDir = '.'; end
    csFile_merge = dir_(vcFile_dir);

    % vsDir = dir(vcFile_dir);
    % vrDatenum = cell2mat_({vsDir.datenum});
    % csFile_merge = {vsDir.name};
    switch fSortByDate
        case 0
        ; %no change
        case 1
        vrDatenum = file_created_meta_(csFile_merge, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);
        case 2
        vrDatenum = file_created_(csFile_merge, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);
        case 3
        vrDatenum = file_created_(csFile_merge, 'lastModifiedTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);
        case 4
        [csFile_merge, ix] = sort_nat_(csFile_merge, 'ascend');
        otherwise
        fprintf(2, 'dir_file_: Invalid option: %d\n', fSortByDate);
    end
    % csFile_merge1 = csFile_merge;
    % csFile_merge = cellfun(@(vc)[vcDir, filesep(), vc], csFile_merge, 'UniformOutput', 0);
end %func
