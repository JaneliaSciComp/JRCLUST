%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function fSuccess = fwrite_(fid, vr)
    try
        fwrite(fid, vr, class(vr));
        fSuccess = 1; % TW
    catch
        fSuccess = 0; % TW
    end
end %func
