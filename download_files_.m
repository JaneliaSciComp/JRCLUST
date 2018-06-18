%--------------------------------------------------------------------------
function vlSuccess = download_files_(csLink, csDest)
    % download file from the web
    if nargin<2, csDest = link2file_(csLink); end
    vlSuccess = false(size(csLink));
    for iFile=1:numel(csLink)
        try
            % download from list of files
            fprintf('\tDownloading %s: ', csLink{iFile});
            vcFile_out1 = websave(csDest{iFile}, csLink{iFile});
            fprintf('saved to %s\n', vcFile_out1);
            vlSuccess(iFile) = 1;
        catch
            fprintf(2, '\n\tCannot download. Check internet connection.\n');
        end
    end %for
end %func
