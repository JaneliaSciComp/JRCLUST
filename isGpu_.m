%--------------------------------------------------------------------------
function flag = isGpu_(vr)
    try
        flag = isa(vr, 'gpuArray');
    catch
        flag = 0;
    end
end
