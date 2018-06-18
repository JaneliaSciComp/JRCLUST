%--------------------------------------------------------------------------
function [vrD1, mrPv1] = cov2var_(mrCov1, nPc)
    [mrPv1, vrD1] = eig(mrCov1);
    vrD1 = cumsum(flipud(diag(vrD1)));
    vrD1 = vrD1 / vrD1(end);
    if nargin>=2
        vrD1 = vrD1(nPc);
    end
end
