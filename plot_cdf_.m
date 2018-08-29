%--------------------------------------------------------------------------
function plot_cdf_(vrSnr, fNorm)
    if nargin<2, fNorm=0; end
    vrX = sort(vrSnr,'ascend');
    vrY = 1:numel(vrSnr);
    if fNorm, vrY=vrY/vrY(end); end
    stairs(vrY, vrX);
end % function
