%--------------------------------------------------------------------------
% 8/7/17 JJJ: tested and documented, file creation date
function [multiFilenames, vcDir] = dir_file_(vcFile_dir, fSortByDate)
    % function [multiFilenames, vcDir, multiFilenames1] = dir_file_(vcFile_dir, fSortByDate)
    % search for files and sort by date

    if nargin<2, fSortByDate=1; end

    [vcDir, ~, ~] = fileparts(vcFile_dir);
    if isempty(vcDir), vcDir = '.'; end
    multiFilenames = dir_(vcFile_dir);

    % vsDir = dir(vcFile_dir);
    % vrDatenum = cell2mat_({vsDir.datenum});
    % multiFilenames = {vsDir.name};
    switch fSortByDate
        case 0
        ; %no change
        case 1
        vrDatenum = file_created_meta_(multiFilenames, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        multiFilenames = multiFilenames(ix);
        case 2
        vrDatenum = file_created_(multiFilenames, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        multiFilenames = multiFilenames(ix);
        case 3
        vrDatenum = file_created_(multiFilenames, 'lastModifiedTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        multiFilenames = multiFilenames(ix);
        case 4
        [multiFilenames, ix] = sort_nat_(multiFilenames, 'ascend');
        otherwise
        fprintf(2, 'dir_file_: Invalid option: %d\n', fSortByDate);
    end
    % multiFilenames1 = multiFilenames;
    % multiFilenames = cellfun(@(vc)[vcDir, filesep(), vc], multiFilenames, 'UniformOutput', 0);
end %func
