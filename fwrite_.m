%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function fSuccess = fwrite_(fid, vr)
    try
        fwrite(fid, vr, class(vr));
        fSucces = 1;
    catch
        fSucces = 0;
    end
end %func
