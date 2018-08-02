%--------------------------------------------------------------------------
function [mrFet1_clu1, iSite_clu1] = S_clu_getFet_(S_clu, iClu, viSite2_spk)
    global trFet_spk
    if nargin<3, viSite2_spk = get0_('viSite2_spk'); end
    iSite_clu1 = S_clu.clusterSites(iClu);
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    if isempty(viSpk_clu1), mrFet1_clu1=[]; return; end
    mrFet1_clu1 = squeeze_(trFet_spk(:,1,viSpk_clu1));
    if size(trFet_spk,2) >= 2
        mrFet2_clu1 = squeeze_(trFet_spk(:,2,viSpk_clu1));
        viSpk_clu1_site2 = find(viSite2_spk(viSpk_clu1) == iSite_clu1);
        mrFet1_clu1(:,viSpk_clu1_site2) = mrFet2_clu1(:,viSpk_clu1_site2);
    end
end %func
