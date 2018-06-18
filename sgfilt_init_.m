%--------------------------------------------------------------------------
function [miA, miB, viC] = sgfilt_init_(nData, nFilt, fGpu)
    persistent miA_ miB_ viC_ nData_ nFilt_
    if nargin<2, fGpu=0; end

    % Build filter coeff
    if isempty(nData_), nData_ = 0; end
    try a = size(miA_); catch, nData_ = 0; end
    if nData_ == nData && nFilt_ == nFilt
        [miA, miB, viC] = deal(miA_, miB_, viC_);
    else
        vi0 = gpuArray_(int32(1):int32(nData), fGpu)';
        vi1 = int32(1):int32(nFilt);
        miA = min(max(bsxfun(@plus, vi0, vi1),1),nData);
        miB = min(max(bsxfun(@plus, vi0, -vi1),1),nData);
        viC = gpuArray_(int32(-nFilt:nFilt), fGpu);
        [nData_, nFilt_, miA_, miB_, viC_] = deal(nData, nFilt, miA, miB, viC);
    end
end %func
