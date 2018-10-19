%--------------------------------------------------------------------------
function [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk)
    % decide whether to use 1, 2, or 3 features
    if nargin<5, vlRedo_spk = []; end
    nFet_use = get_set_(P, 'nFet_use', 2);

    time_feature_factor = get_set_(P, 'time_feature_factor', 0); % TW

    [viSpk1_, viSpk2_] = multifun_(@int32, S0.cviSpk_site{iSite}, S0.cviSpk2_site{iSite});
    if ~isempty(vlRedo_spk)
        viSpk1_ = viSpk1_(vlRedo_spk(viSpk1_));
        viSpk2_ = viSpk2_(vlRedo_spk(viSpk2_));
    end
    if isempty(viSpk1_)
        [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = deal([]);
        return;
    end

    if isempty(viSpk2_), nFet_use=1; end

    switch nFet_use
        case 3
        viSpk3_ = int32(S0.cviSpk3_site{iSite});
        if ~isempty(viSpk3_), viSpk3_ = viSpk3_(vlRedo_spk(viSpk3_)); end
        mrFet12_ = [squeeze_(trFet_spk(:,1,viSpk1_),2), squeeze_(trFet_spk(:,2,viSpk2_),2), squeeze_(trFet_spk(:,3,viSpk3_),2)];
        viSpk12_ = [viSpk1_; viSpk2_; viSpk3_];
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_) + numel(viSpk3_));
        case 2
        mrFet12_ = [squeeze_(trFet_spk(:,1,viSpk1_),2), squeeze_(trFet_spk(:,2,viSpk2_),2); single(S0.viTime_spk([viSpk1_;viSpk2_]))']; % TW
        mrFet12_(end, :) = time_feature_factor*std(mrFet12_(1, :)).*mrFet12_(end, :)./std(mrFet12_(end, :)); % TW
        viSpk12_ = [viSpk1_; viSpk2_];
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
        case 1
        mrFet12_ = [squeeze_(trFet_spk(:,1,viSpk1_),2); single(S0.viTime_spk(viSpk1_))']; % TW
        mrFet12_(end, :) = time_feature_factor*std(mrFet12_(1, :)).*mrFet12_(end, :)./std(mrFet12_(end, :)); % TW
        viSpk12_ = viSpk1_;
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
    end

    try
        nSites_fet = 1 + P.maxSite*2 - P.nSites_ref;
        if get_set_(P, 'fSpatialMask_clu', 1) && nSites_fet >= get_set_(P, 'min_sites_mask', 5)
            nFetPerChan = size(mrFet12_,1) / nSites_fet;
            vrSpatialMask = spatialMask_(P, iSite, nSites_fet, P.maxDist_site_um);
            vrSpatialMask = repmat(vrSpatialMask(:), [nFetPerChan, 1]);
            mrFet12_ = bsxfun(@times, mrFet12_, vrSpatialMask(:));
        end
    catch
        disperr_('Spatial mask error');
    end
    if get_set_(P, 'fSqrt_fet', 0), mrFet12_ = signsqrt_(mrFet12_); end
    if get_set_(P, 'fLog_fet', 0), mrFet12_ = signlog_(mrFet12_); end
    if get_set_(P, 'fSquare_fet', 0), mrFet12_ = (mrFet12_).^2; end
    viiSpk12_ord_ = rankorder_(viSpk12_, 'ascend');
end %func
