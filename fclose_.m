%--------------------------------------------------------------------------
function fid = fclose_(fid, fVerbose)
    % Sets fid = [] if closed properly
    if nargin < 2
        fVerbose = 0;
    end

    if isempty(fid) || ischar(fid)
        return;
    end

    try
        fclose(fid);
        fid = [];
        if fVerbose
            disp('File closed.');
        end
    catch
        disperr_();
    end
end %func
