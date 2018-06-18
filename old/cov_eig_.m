%--------------------------------------------------------------------------
function dist12 = cov_eig_(mrCov1, mrCov2, nPc)
    d1 = cov2var_(mrCov1, nPc);
    d2 = cov2var_(mrCov2, nPc);
    d12 = cov2var_(mrCov1 + mrCov2, nPc);
    dist12 = d12 / ((d1+d2)/2);
end %func
