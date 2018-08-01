%--------------------------------------------------------------------------
% 8/11/17 JJJ: Initial implementation
function S_fig = Fig_preview_update_(hFig, S_fig, fKeepView)
    % Handle parameter update requests
    % S_fig = Fig_preview_update_(hFig, S_fig, P)
    % S_fig = Fig_preview_update_(hFig, S_fig, fUpdate)

    if nargin<1, hFig = getCachedFig('Fig_preview'); end
    if nargin<2, S_fig = get(hFig, 'UserData'); end
    if nargin<3, fKeepView = 0; end

    P = get0_('P');
    nSites = numel(P.chanMap);
    figure_wait_(1, hFig); drawnow;
    fft_thresh = S_fig.fft_thresh;
    if fft_thresh > 0
        S_fig.mnWav_clean = fft_clean_(S_fig.mnWav_raw, struct_add_(P, fft_thresh)); % fft filter
    else
        S_fig.mnWav_clean = S_fig.mnWav_raw;
    end

    % Find bad sites
    if S_fig.thresh_corr_bad_site > 0
        S_fig.vlSite_bad = S_fig.vrCorr_max_site < S_fig.thresh_corr_bad_site;
        S_fig.viSite_bad = find(S_fig.vlSite_bad);
    elseif ~isempty(S_fig.viSiteZero)
        S_fig.viSite_bad = S_fig.viSiteZero;
        S_fig.vlSite_bad = false(size(S_fig.vrCorr_max_site));
        S_fig.vlSite_bad(S_fig.viSiteZero) = 1;
        S_fig.viSiteZero = [];
    else
        S_fig.vlSite_bad = false(size(S_fig.vrCorr_max_site));
        S_fig.viSite_bad = [];
    end

    % Perform filter fft_thresh
    P_ = set_(P, 'vcCommonRef', 'none', 'useGPU', 0, 'vcFilter', S_fig.vcFilter, ...
    'blank_period_ms', S_fig.blank_period_ms, 'blank_thresh', S_fig.blank_thresh, 'fParfor', 0);
    mnWav_filt = filt_car_(S_fig.mnWav_clean, P_);
    % if strcmpi(S_fig.vcCommonRef, 'median')
    %     vrWav_filt_mean = median(mnWav_filt(:,~S_fig.vlSite_bad), 2);
    % else
    %     vrWav_filt_mean = mean(mnWav_filt(:,~S_fig.vlSite_bad), 2);
    % end
    vrWav_filt_mean = mr2ref_(mnWav_filt, S_fig.vcCommonRef, S_fig.viSite_bad);
    if ~strcmpi(S_fig.vcCommonRef, 'none')
        mnWav_filt = bsxfun(@minus, mnWav_filt, int16(vrWav_filt_mean));
    end
    S_fig.vrWav_filt_mean = madscore_(mean(vrWav_filt_mean, 2)); % Save in MAD unit
    S_fig.mnWav_filt = mnWav_filt;
    % reference threshold. determine
    % vrWav_filt_mean = mean(mnWav_filt, 2) * P.uV_per_bit; % @TODO: save in MAD unit

    % vrPower_psd = abs(mean(fft(S_fig.mnWav_raw(:,~S_fig.vlSite_bad)), 2));
    [mrPower_psd, S_fig.vrFreq_psd] = psd_(S_fig.mnWav_raw(:,~S_fig.vlSite_bad), P.sRateHz, 4);
    [mrPower_clean_psd] = psd_(S_fig.mnWav_clean(:,~S_fig.vlSite_bad), P.sRateHz, 4);
    [S_fig.vrPower_psd, S_fig.vrPower_clean_psd] = multifun_(@(x)mean(x,2), mrPower_psd, mrPower_clean_psd);

    % Apply threshold and perform spike detection
    vrRmsQ_site = mr2rms_(mnWav_filt, 1e5);
    vnThresh_site = int16(vrRmsQ_site * S_fig.qqFactor);
    vnThresh_site(S_fig.vlSite_bad) = nan; % shows up as 0 for int16
    S_fig.mlWav_thresh = bsxfun(@lt, mnWav_filt, -abs(vnThresh_site)); %negative threshold crossing
    S_fig.mlWav_thresh(:, S_fig.vlSite_bad) = 0;
    S_fig.vnThresh_site = vnThresh_site;

    % Spike detection
    % P_.fMerge_spk = 0;
    [vlKeep_ref, S_fig.vrMad_ref] = car_reject_(vrWav_filt_mean, P_);
    [S_fig.spikeTimes, S_fig.vnAmp_spk, viSite_spk] = detect_spikes_(mnWav_filt, vnThresh_site, vlKeep_ref, P_);
    t_dur = size(mnWav_filt,1) / P.sRateHz;
    S_fig.vrEventRate_site = hist(viSite_spk, 1:nSites) / t_dur; % event count
    S_fig.vrEventSnr_site = abs(single(arrayfun(@(i)median(S_fig.vnAmp_spk(viSite_spk==i)), 1:nSites))) ./ vrRmsQ_site;
    S_fig.vlKeep_ref = vlKeep_ref;
    S_fig.viSite_spk = viSite_spk;
    % Spike stats: such as # sites/event over threshold

    % Exit
    set(hFig, 'UserData', S_fig);
    [hFig, S_fig] = Fig_preview_plot_(P, fKeepView);
end %func
