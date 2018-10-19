%--------------------------------------------------------------------------
function [mr, fGpu] = gpuArray_(mr, fGpu)
    if nargin<2, fGpu = 1; end
    if ~fGpu, return; end
    try
        mr = gpuArray(mr);
        fGpu = 1;
    catch
        try % retry after resetting the GPU memory
            gpuDevice(1);
            mr = gpuArray(mr);
            fGpu = 1;
        catch % no GPU device found
            fGpu = 0;
        end
    end
end
