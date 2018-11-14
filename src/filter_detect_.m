%--------------------------------------------------------------------------
function [mn1, nShift_post] = filter_detect_(mn, P, vcMode)
    % returns spatial sd
    % mn0 = single(mn);
    % mn0 = bsxfun(@minus, mn0, mean(mn0, 2)) .^ 2;
    % 6/29/17 JJJ: filter detection

    if nargin<3, vcMode = get_(P, 'vcFilter_detect'); end

    viSites_use = 1:(1+2*P.maxSite - P.nSites_ref);
    viSites_ref = (1+2*P.maxSite - P.nSites_ref+1):(1+2*P.maxSite);
    fprintf('filter_detect\n\t'); t1= tic;
    miSites = jrclust.utils.tryGpuArray(P.miSites(viSites_use, :));
    miSites_ref = jrclust.utils.tryGpuArray(P.miSites(viSites_ref, :));
    nShift_post = 0;
    switch lower(vcMode)
        case 'ndist'
        mn1 = ndist_filt_(mn, get_set_(P, 'ndist_filt', 5));
        case 'chancor'
        mn1 = chancor_(mn, P);
        case 'matched'
        vrFilt_spk = get0_('vrFilt_spk');
        if isempty(vrFilt_spk)
            lim_ = round([3,5]/8 * size(mn,1));
            mn_ = mn(lim_(1):lim_(2),:);
            [vrFilt_spk, vrVaf, nShift_post] = calc_matched_filt_(mn_, P); %detect primary
            set0_(vrFilt_spk, nShift_post);
        else
            nShift_post = get0_('nShift_post');
        end
        %         nShift_post = nShift_post_;
        mn1 = int16(conv2(single(jrclust.utils.tryGather(mn)), vrFilt_spk(:), 'same'));
        %         mn1 = shift_mr_(mn1, nShift_post); % or do it later in the spike detection phase
        %         figure; plot(xcorr(mn(:,41), mn1(:,41), 10));
        %         nShift_post = P.spkLim(1)-1;
        case 'autocov'
        mn1 = -int16(filt_corr(single(mn), 2));
        case 'std-chan'
        mn1 = zeros(size(mn), 'like', mn);
        for iSite=1:size(mn,2)
            %mn_ = mn(:, P.miSites(viSites_use, iSite));
            %vn_ref = mean(mn_(:,viSites_ref),2);
            %mn_ = bsxfun(@minus, mn_(:,viSites_use), vn_ref);
            % mn1(:, iSite) = -int16(std(single(mn(:, P.miSites(:, iSite))), 0, 2));
            mn1(:, iSite) = -int16(std(single(mn(:, miSites(:,iSite))), 1, 2));
            %mn1(:,iSite) = mean(mn_.^2,2) - mean(mn_,2).^2;
            fprintf('.');
        end

        case 'std-time'
        % envelop filter. this affects the threshold.

        case 'nmean'
        mn1 = zeros(size(mn), 'like', mn);
        for iSite=1:size(mn,2)
            %mn_ = mn(:, P.miSites(viSites_use, iSite));
            %vn_ref = mean(mn_(:,viSites_ref),2);
            %mn_ = bsxfun(@minus, mn_(:,viSites_use), vn_ref);
            % mn1(:, iSite) = -int16(std(single(mn(:, P.miSites(:, iSite))), 0, 2));
            %mn1(:, iSite) = mn(:,iSite) - int16(mean(single(mn(:,miSites_ref(:,iSite))),2));
            mn1(:, iSite) = int16(mean(mn(:,miSites(:,iSite)), 2));
            %mn1(:,iSite) = mean(mn_.^2,2) - mean(mn_,2).^2;
            fprintf('.');
        end

    end %switch
    fprintf('\n\ttook %0.1fs\n', toc(t1));
    % mn1 = -int16(sqrt(mn1 / size(P.miSites,1)));
end %func
