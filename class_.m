%--------------------------------------------------------------------------
function [vc, useGPU] = class_(vr)
    % Return the class for GPU or CPU arrays
    if isempty(vr)
        vc = class(gather_(vr));
    else
        vc = class(gather_(vr(1)));
    end
    if nargout>=2, useGPU = isGpu_(vr); end
end % function
