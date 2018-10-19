%--------------------------------------------------------------------------
function dc2_ = compute_dc2_(mrFet12, viiSpk12_ord, n1_, n2_, P)
    % subsample
    % fSubsample12 = 1;
    if get_set_(P, 'fDc_spk', 0)
        dc2_ = (P.dc_percent/100).^2;
        return;
    end % spike-specific dc
    if 0 % 122717 JJJ
        [n1_max, n2_max] = deal(P.dc_subsample, P.dc_subsample * 40);
    else % old
        [n1_max, n2_max] = deal(P.dc_subsample, P.dc_subsample * 4);
    end
    vi1_ = subsample_vr_(1:n1_, n1_max);
    viiSpk1_ord_ = viiSpk12_ord(vi1_);
    mrFet1_ = mrFet12(:,vi1_);

    % if fSubsample12
    %     vi12_ = subsample_vr_(1:(n1_+n2_), P.dc_subsample * 4);
    vi12_ = subsample_vr_(1:n1_, n2_max);
    viiSpk12_ord = viiSpk12_ord(vi12_);
    mrFet12 = mrFet12(:,vi12_);
    % end
    for iRetry=1:2
        try
            mrDist11_2_ = eucl2_dist_(mrFet12, mrFet1_);
            % if get_set_(P, 'f_dpclus', 1)
            mlKeep11_2_ = abs(bsxfun(@minus, viiSpk12_ord, viiSpk1_ord_')) < (n1_+n2_) / P.nTime_clu;
            mrDist11_2_(~mlKeep11_2_) = nan;
        catch
            mrFet12 = gather_(mrFet12);
        end
    end
    % end
    if get_set_(P, 'fDc_subsample_mode', 0)
        mrDist11_2_(mrDist11_2_<=0) = nan;
        dc2_ = quantile(mrDist11_2_(~isnan(mrDist11_2_)), P.dc_percent/100);
    else
        %     mrDist_sub = subsample_mr_(mrDist11_2_, P.dc_subsample, 2);
        mrDist_sub = gather_(mrDist11_2_);
        mrDist_sub(mrDist_sub<=0) = nan;
        if 1
            dc2_ = nanmedian(quantile(mrDist_sub, P.dc_percent/100));
        else
            dc2_ = nan;
        end
        if isnan(dc2_), dc2_ = quantile(mrDist_sub(:), P.dc_percent/100); end
    end
end %func
