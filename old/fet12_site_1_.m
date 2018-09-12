%--------------------------------------------------------------------------
function [mrFet1_, mrFet2_, viSpk1_, viSpk2_, n1_, n2_, viiSpk12_ord_] = fet12_site_1_(mrFet1, mrFet2, S0, P, iSite)

    [viSpk1_, viSpk2_] = deal(S0.cviSpk_site{iSite}, S0.cviSpk2_site{iSite});
    [mrFet1_, mrFet2_] = deal(mrFet1(:,viSpk1_), mrFet2(:,viSpk2_));
    [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
    if P.fGpu
        [mrFet1_, mrFet2_, viSpk1_, viSpk2_] = multifun_(@gpuArray_, mrFet1_, mrFet2_, viSpk1_, viSpk2_);
    end
    viiSpk12_ord_ = rankorder_([viSpk1_; viSpk2_], 'ascend');
end %func
