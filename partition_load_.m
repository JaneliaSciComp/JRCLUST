%--------------------------------------------------------------------------
function [nLoad1, nSamples_load1, nSamples_last1] = partition_load_(nSamples1, nSamples_max)
    nLoad1 = setlim_(ceil(nSamples1 / nSamples_max), [1, inf]);
    nSamples_load1 = min(nSamples1, nSamples_max);
    if nLoad1 == 1
        nSamples_load1 = nSamples1;
        nSamples_last1 = nSamples1;
    else
        nSamples_last1 = mod(nSamples1, nSamples_load1);
        if nSamples_last1==0
            nSamples_last1 = nSamples_load1;
        elseif nSamples_last1 < nSamples_load1/2
            % if last part is too small increase the size
            nLoad1 = nLoad1 - 1;
            nSamples_last1 = nSamples_last1 + nSamples_load1;
        end
    end
end % function
