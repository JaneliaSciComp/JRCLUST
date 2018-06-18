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
