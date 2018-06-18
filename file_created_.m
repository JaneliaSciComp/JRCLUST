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
