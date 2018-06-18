%--------------------------------------------------------------------------
% 9/29/17 JJJ: delta==0 (found other spike that has exactly same PCA score) is set to nan
% 9/5/17 JJJ: Created
function [icl, x, z] = detrend_local_(S_clu, P, fLocal)
    if nargin<3, fLocal = 1; end
    maxCluPerSite = get_set_(P, 'maxCluPerSite', 20); % get 10 clu per site max
    S0 = get0();
    % cvi_rho_rank_site = cellfun(@(vi)rankorder_(S_clu.rho(vi), 'ascend'), S0.cviSpk_site, 'UniformOutput', 0);
    % cvi_delta_rank_site = cellfun(@(vi)rankorder_(S_clu.delta(vi), 'ascend'), S0.cviSpk_site, 'UniformOutput', 0);

    % detrend for each site and apply
    if fLocal
        cvi_rho_site = cellfun(@(vi)S_clu.rho(vi), S0.cviSpk_site, 'UniformOutput', 0);
        cvi_delta_site = cellfun(@(vi)S_clu.delta(vi), S0.cviSpk_site, 'UniformOutput', 0);
        %     [x0, y0] = deal(S_clu.rho, S_clu.delta);
        cvi_cl = cell(size(S0.cviSpk_site));
        x = log10(S_clu.rho);
        z = zeros(size(S_clu.delta), 'like', S_clu.delta);
        for iSite = 1:numel(S0.cviSpk_site)
            [rho_, delta_] = deal(cvi_rho_site{iSite}, cvi_delta_site{iSite});
            x_ = log10(rho_);
            %         y_ = log10(rho_) + log10(delta_);
            y_ = (delta_);
            viDetrend = find(delta_ < 1 & delta_ > 0 & rho_ > 10^P.rho_cut & rho_ < .1 & isfinite(x_) & isfinite(y_));
            [y_det_, z_] = detrend_(x_, y_, viDetrend, 1);
            viSpk_ = S0.cviSpk_site{iSite};
            if isempty(viSpk_), continue; end
            vl_zero_ = find(delta_==0);
            [y_det_(vl_zero_), z_(vl_zero_)] = nan;
            [icl_, vrZ_] = find_topn_(y_det_, maxCluPerSite, ...
            find(rho_ > 10^P.rho_cut & ~isnan(y_det_)));
            if isempty(icl_), continue; end
            cvi_cl{iSite} = viSpk_(icl_);
            z(viSpk_) = z_;
            %         figure; plot(x_, y_det_, 'k.', x_(icl_), y_det_(icl_), 'ro'); axis([-5 0 -1 1]); grid on;
        end
        icl = cell2vec_(cvi_cl);
    else
        x = log10(S_clu.rho);
        y = S_clu.delta;
        viDetrend = find(S_clu.delta < 1 & S_clu.delta > 0 & S_clu.rho > 10^P.rho_cut & S_clu.rho < .1 & isfinite(x) & isfinite(y));
        [~, z] = detrend_(x, y, viDetrend, 1);
        z(S_clu.delta==0) = nan;
        [icl, vrZ_] = find_topn_(z, maxCluPerSite * numel(S0.cviSpk_site), ...
        find(S_clu.rho > 10^P.rho_cut & ~isnan(z)));
        icl(vrZ_ < 10^P.delta1_cut) = [];
        %     icl = find(x>=P.rho_cut & z>=10^P.delta1_cut);
    end


    if nargout==0
        figure; plot(x,z,'.', x(icl),z(icl),'ro'); grid on;
        axis_([-5 0 -20 100]);
        title(sprintf('%d clu', numel(icl)));
    end
end %func
