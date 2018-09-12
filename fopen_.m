%--------------------------------------------------------------------------
function [fid, nBytes, header_offset] = fopen_(vcFile, vcMode)
    if nargin < 2, vcMode = 'r'; end
    try
        if matchFileExt_(vcFile, {'.ns5', '.ns2'})
            [fid, nBytes, header_offset] = fopen_nsx_(vcFile);
        else
            fid = fopen(vcFile, vcMode);
            nBytes = getBytes_(vcFile);
            header_offset = 0;
        end
    catch
        disperr_();
        fid = [];
        nBytes = [];
    end
end %func
