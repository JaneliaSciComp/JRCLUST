%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function fSuccess = fwrite_(fid, vr)
    try
        fwrite(fid, vr, class(vr));
        fSuccess = 1;
    catch
        fSuccess = 0;
    end
end %func
