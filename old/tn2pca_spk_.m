%--------------------------------------------------------------------------
function mrPv = tn2pca_spk_(tn, nPc)
    % calc cov by averaging and compute pca
    [nT, nC, nSpk] = size(tn);
    % mrCov1 = zeros(nT, nT, 'single');
    mr = single(reshape(tn, nT,[]));
    mr = meanSubt_(mr);
    [~, mrPv] = pca(mr, 'NumComponents', nPc);
    % tr = meanSubt_(single(tn));
    % for iSpk = 1:nSpk
    %     mr1 = tr(:,:,iSpk);
    %     mrCov1 = mrCov1 + mr1 * mr1';
    % end
    % [mrPv, vrL] = eig(mrCov1);
    % mrPv = mrPv(:, end:-1:end-nPc+1);
end %func
