%--------------------------------------------------------------------------
% 8/7/17 JJJ: tested and documented, file creation date
function [csFile_merge, vcDir] = dir_file_(vcFile_dir, fSortByDate)
    % function [csFile_merge, vcDir, csFile_merge1] = dir_file_(vcFile_dir, fSortByDate)
    % search for files and sort by date

    if nargin < 2
        fSortByDate=1;
    end

    [vcDir, ~, ~] = fileparts(vcFile_dir);
    if isempty(vcDir)
        vcDir = '.';
    end

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
            [csFile_merge, ix] = jrclust.utils.sortNat(csFile_merge, 'ascend');
        otherwise
            fprintf(2, 'dir_file_: Invalid option: %d\n', fSortByDate);
    end
    % csFile_merge1 = csFile_merge;
    % csFile_merge = cellfun(@(vc)[vcDir, filesep(), vc], csFile_merge, 'UniformOutput', 0);
end %func

%% local functions

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
end %func


%--------------------------------------------------------------------------
% 8/7/17 JJJ: Created, tested and documented
function vrDatenum = file_created_(csFiles, vcDateMode)
    % Return the file creation time. Matlab "dir" command returns last modified time
    % vcDateMode: {'lastModifiedTime', 'creationTime', 'lastAccessTime'};
    if nargin<2, vcDateMode = 'creationTime'; end

    if ischar(csFiles), csFiles = {csFiles}; end
    vrDatenum = nan(size(csFiles));
    for iFile=1:numel(csFiles)
        vcFile_ = csFiles{iFile};
        try
            vcDate = char(java.nio.file.Files.readAttributes(java.io.File(vcFile_).toPath(), vcDateMode, javaArray('java.nio.file.LinkOption', 0)));
            iStart = find(vcDate=='=', 1, 'first');
            iEnd = find(vcDate=='Z', 1, 'last');
            vcDate = vcDate(iStart+1:iEnd-1);
            vcDate(vcDate=='T') = ' '; % replace with blank
            vrDatenum(iFile) = datenum(vcDate, 'yyyy-mm-dd HH:MM:SS.FFF');
        catch
            fprintf(2, 'File not found: %d\n', vcFile_);
        end
    end
    % datestr(datenum1, 'yyyy-mm-dd HH:MM:SS.FFF')
end %func

%--------------------------------------------------------------------------
% 8/7/17 JJJ: Created, tested and documented
function vrDatenum = file_created_meta_(csFiles, vcDateMode)
    % Read the SpikeGLX meta file for the file cration time if available
    if nargin<2, vcDateMode = 'creationTime'; end
    if ischar(csFiles), csFiles = {csFiles}; end
    vrDatenum = nan(size(csFiles));
    for iFile=1:numel(csFiles)
        vcFile_ = csFiles{iFile};
        try
            vcFile_meta_ = subsFileExt_(vcFile_, '.meta');
            if exist(vcFile_meta_, 'file')
                S_meta_ = text2struct_(vcFile_meta_);
                vcDatenum_ = S_meta_.fileCreateTime;
                vcDatenum_(vcDatenum_=='T') = ' ';
                vrDatenum(iFile) = datenum(vcDatenum_, 'yyyy-mm-dd HH:MM:SS');
            end
        catch
            ;
        end
        if isnan(vrDatenum(iFile))
            vrDatenum(iFile) = file_created_(vcFile_, vcDateMode);
        end
    end %for
end %func
