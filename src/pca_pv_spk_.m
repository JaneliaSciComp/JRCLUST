%--------------------------------------------------------------------------
function [mrPv1, mrPv2] = pca_pv_spk_(sampledSpikes, iSite, sampledWindows)
    % if viSite not found just set to zero
    nSites = numel(iSite);
    S0 = get0_();
    mrPv_global = get_(S0, 'mrPv_global');
    if ~isempty(mrPv_global)
        % show global pca
        mrPv1 = repmat(mrPv_global(:,1), [1, nSites]);
        mrPv2 = repmat(mrPv_global(:,2), [1, nSites]);
    else
        if nargin<3
            sampledWindows = permute(jrclust.utils.getSampledWindows(sampledSpikes, iSite, S0, 0), [1, 3, 2]);
        end
        nT = size(sampledWindows, 1);
        % show site pca
        [mrPv1, mrPv2] = deal(zeros(nT, nSites, 'single'));
        for jSite = 1:nSites
            [mrPv1(:, jSite), mrPv2(:, jSite)] = pca_pv_(sampledWindows(:, :, jSite));
        end %for
    end
end %func
