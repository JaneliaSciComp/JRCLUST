%--------------------------------------------------------------------------
function nBytes = mem_max_(P)
    try
        if P.fGpu
            S = gpuDevice(); % does not reset GPU
            nBytes = floor(S(1).AvailableMemory());
        else
            S = memory();
            nBytes = floor(S.MaxPossibleArrayBytes());
        end
    catch
        nBytes = inf; % assume infinite memory
    end
end %func
