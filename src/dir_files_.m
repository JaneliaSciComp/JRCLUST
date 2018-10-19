%--------------------------------------------------------------------------
function csFiles_bin = dir_files_(csFile_merge, vcFile_txt, file_sort_merge)
    % Display binary files
    if nargin<2, vcFile_txt=''; end
    if nargin<3, file_sort_merge = []; end
    if isempty(file_sort_merge)
        file_sort_merge = 1;
    elseif ischar(file_sort_merge)
        file_sort_merge = str2double(file_sort_merge);
    end

    if isempty(csFile_merge)
        fprintf('No files to merge ("csFile_merge" is empty).\n');
        csFiles_bin = {}; return;
    end
    % fprintf('Listing files to merge ("csFile_merge"):\n');
    csFiles_bin = filter_files_(csFile_merge, file_sort_merge);
    if nargout==0
        %     arrayfun(@(i)fprintf('%d: %s\n', i, csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
        arrayfun(@(i)fprintf('%s\n', csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
    end
    if ~isempty(vcFile_txt)
        cellstr2file_(vcFile_txt, csFiles_bin);
        fprintf('%s is created.\n', vcFile_txt);
        edit(vcFile_txt);
    end
end %func
