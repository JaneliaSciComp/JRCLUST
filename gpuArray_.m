%--------------------------------------------------------------------------
function [mr, useGPU] = gpuArray_(mr, useGPU)
    if nargin<2, useGPU = 1; end
    if ~useGPU, return; end
    try
        mr = gpuArray(mr);
        useGPU = 1;
    catch
        try % retry after resetting the GPU memory
            gpuDevice(1);
            mr = gpuArray(mr);
            useGPU = 1;
        catch % no GPU device found
            useGPU = 0;
        end
    end
end
