function [val, success] = tryGpuArray(val, useGPU)
    %TRYGPUARRAY Try to construct a gpuArray from val
    if nargin == 2 && ~useGPU
        return;
    end

    try
        val = gpuArray(val);
        success = 1;
    catch
        try % reset GPU memory and try again
            gpuDevice(1);
            val = gpuArray(val);
            success = 1;
        catch % no GPU device found
            success = 0;
        end
    end
end
