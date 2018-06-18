%--------------------------------------------------------------------------
function mrPv = tn2pca_cov_(tn, nPc)
    % calc cov by averaging and compute pca
    [nT, nC, nSpk] = size(tn);
    mrCov1 = zeros(nT, nT, 'double');
    tr = meanSubt_(single(tn));
    for iSpk = 1:nSpk
        mr1 = tr(:,:,iSpk);
        mrCov1 = mrCov1 + double(mr1 * mr1');
    end
    [mrPv, vrL] = eig(mrCov1);
    mrPv = mrPv(:, end:-1:end-nPc+1);
end %func
