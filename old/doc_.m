%--------------------------------------------------------------------------
% 10/8/17 JJJ: Opens the doc file from the current JRC folder
% 7/24/17 JJJ: open the latest manual from Dropbox if the link is valid
function doc_(vcFile_doc)
    % Open JRCLUST PDF documentation
    if nargin<1, vcFile_doc = 'JRCLUST manual.pdf'; end
    vcFile_doc = jrcpath_(vcFile_doc);
    if exist_file_(vcFile_doc)
        disp(vcFile_doc);
        open(vcFile_doc);
    else
        fprintf(2, 'File does not exist: %s\n', vcFile_doc);
    end
end %func
