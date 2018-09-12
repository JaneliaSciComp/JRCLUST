%--------------------------------------------------------------------------
function [mrPv1, mrPv2] = pca_pv_spk_(viSpk1, viSites1, tnWav_spk1)
    % if viSite not found just set to zero
    nSites = numel(viSites1);
    S0 = get0_();
    mrPv_global = get_(S0, 'mrPv_global');
    if ~isempty(mrPv_global)
        % show global pca
        mrPv1 = repmat(mrPv_global(:,1), [1, nSites]);
        mrPv2 = repmat(mrPv_global(:,2), [1, nSites]);
    else
        if nargin<3
            tnWav_spk1 = permute(tnWav_spk_sites_(viSpk1, viSites1, S0, 0), [1,3,2]);
        end
        nT = size(tnWav_spk1, 1);
        % show site pca
        [mrPv1, mrPv2] = deal(zeros(nT, nSites, 'single'));
        for iSite1=1:nSites
            [mrPv1(:,iSite1), mrPv2(:,iSite1)] = pca_pv_(tnWav_spk1(:,:,iSite1));
        end %for
    end
end %func
