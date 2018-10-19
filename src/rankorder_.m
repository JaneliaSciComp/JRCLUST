%--------------------------------------------------------------------------
function [vi, viSort] = rankorder_(vr, vcOrder)
    % warning: 32 bit addressing
    if nargin<2, vcOrder = 'ascend'; end
    n=numel(vr);
    [~,viSort] = sort(vr, vcOrder);
    if isGpu_(vr)
        vi = zeros(n,1,'int32', 'gpuArray');
        vi(viSort) = 1:n;
    else
        vi=zeros(n,1,'int32');
        vi(viSort) = 1:n;
    end
end %func
