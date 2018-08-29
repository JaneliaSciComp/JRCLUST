%--------------------------------------------------------------------------
function [mrPatchX, mrPatchY] = probe_map_(P)
    vrX = [0 0 1 1] * P.vrSiteHW(2);
    vrY = [0 1 1 0] * P.vrSiteHW(1);
    mrPatchX = bsxfun(@plus, P.mrSiteXY(:,1)', vrX(:));
    mrPatchY = bsxfun(@plus, P.mrSiteXY(:,2)', vrY(:));
end % function
