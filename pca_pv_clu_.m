%--------------------------------------------------------------------------
function [mrPv1, mrPv2] = pca_pv_clu_(viSites, iClu1, iClu2)
    % [mrPv1, mrPv2] = pca_pv_clu_(viSites, iClu1): return two prinvec per site from same clu
    % [mrPv1, mrPv2] = pca_pv_clu_(viSites, iClu1, iClu2): return one prinvec per clu per site

    if nargin<3, iClu2 = []; end
    S0 = get0_();
    S_clu = S0.S_clu;
    % show site pca
    MAX_SAMPLE = 10000; %for pca
    viSpk1 = subsample_vr_(S_clu.cviSpk_clu{iClu1}, MAX_SAMPLE);
    spikeWaveforms1 = permute(spikeWaveforms_sites_(viSpk1, viSites, S0, 0), [1,3,2]);
    [nT, nSites] = deal(size(spikeWaveforms1, 1), numel(viSites));
    [mrPv1, mrPv2] = deal(zeros(nT, nSites, 'single'));

    if isempty(iClu2)
        for iSite1=1:nSites
            [mrPv1(:,iSite1), mrPv2(:,iSite1)] = pca_pv_(spikeWaveforms1(:,:,iSite1));
        end %for
    else
        viSpk2 = subsample_vr_(S_clu.cviSpk_clu{iClu2}, MAX_SAMPLE);
        spikeWaveforms2 = permute(spikeWaveforms_sites_(S_clu.cviSpk_clu{iClu2}, viSites, S0, 0), [1,3,2]);
        for iSite1=1:nSites
            mrPv1(:,iSite1) = pca_pv_(spikeWaveforms1(:,:,iSite1));
            mrPv2(:,iSite1) = pca_pv_(spikeWaveforms2(:,:,iSite1));
        end %for
    end
end %func
