%--------------------------------------------------------------------------
function [mrFet1_clu1, iSite_clu1] = S_clu_getFet_(S_clu, iClu, spikeSecondarySites)
    global spikeFeatures

    S0 = get(0, 'UserData');

    if ~allDimsEqual(spikeFeatures, S0.featureDims)
        spikeFeatures = getSpikeFeatures();
    end

    if nargin<3, spikeSecondarySites = get0_('spikeSecondarySites'); end
    iSite_clu1 = S_clu.clusterSites(iClu);
    viSpk_clu1 = S_clu.spikesByCluster{iClu};
    if isempty(viSpk_clu1), mrFet1_clu1=[]; return; end
    mrFet1_clu1 = squeeze_(spikeFeatures(:,1,viSpk_clu1));
    if size(spikeFeatures,2) >= 2
        mrFet2_clu1 = squeeze_(spikeFeatures(:,2,viSpk_clu1));
        viSpk_clu1_site2 = find(spikeSecondarySites(viSpk_clu1) == iSite_clu1);
        mrFet1_clu1(:,viSpk_clu1_site2) = mrFet2_clu1(:,viSpk_clu1_site2);
    end
end % function
