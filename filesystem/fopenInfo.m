%--------------------------------------------------------------------------
function [fid, nBytes, headerOffset] = fopenInfo(filename, permission)
    if nargin < 2
        permission = 'r';
    end

    try
        if matchFileExt_(filename, {'.ns5', '.ns2'})
            [fid, nBytes, headerOffset] = fopenNeuroshare(filename);
        else
            fid = fopen(filename, permission);
            nBytes = getBytes_(filename);
            headerOffset = 0;
        end
    catch
        disperr_();
        fid = [];
        nBytes = [];
    end
end %func
