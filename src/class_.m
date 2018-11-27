%--------------------------------------------------------------------------
function [vc, fGpu] = class_(vr)
    % Return the class for GPU or CPU arrays
    if isempty(vr)
        vc = class(jrclust.utils.tryGather(vr));
    else
        vc = class(jrclust.utils.tryGather(vr(1)));
    end
    if nargout >=2
        fGpu = isa(vr, 'gpuArray');
    end
end %func
