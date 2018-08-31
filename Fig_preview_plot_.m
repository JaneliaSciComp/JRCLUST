%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
% @TODO: run spike detection and show detected spikes, or below threshold values
function [hFig, S_fig] = Fig_preview_plot_(P, fKeepView)

    if nargin<1, P = []; end
    if nargin<2, fKeepView = 0; end
    hWait = msgbox_('Plotting...', 0, 1);
    if isempty(P), P = get0_('P'); end
    [hFig, S_fig] = getCachedFig('Fig_preview');

    figure_wait_(1, hFig);
    nSites = size(S_fig.mnWav_filt,2);
    viPlot = S_fig.nlim_bin(1):S_fig.nlim_bin(2);
    vrTime_sec = viPlot / P.sampleRateHz;
    tlim_sec = (S_fig.nlim_bin + [-1 1]) / P.sampleRateHz;


    %-----
    % Mean plot
    if isempty(get_(S_fig, 'hPlot_mean'))
        S_fig.hPlot_mean = plot(S_fig.hAx_mean, nan, nan, 'k');
        S_fig.hPlot_mean_thresh = plot(S_fig.hAx_mean, nan, nan, 'r');
    end
    switch S_fig.vcRef_view
        case 'original'
        vrY_ = S_fig.vrWav_filt_mean(viPlot);
        case 'binned'
        vrY_ = S_fig.vrMad_ref(viPlot);
        otherwise, disperr_();
    end
    set(S_fig.hPlot_mean, 'XData', vrTime_sec, 'YData', vrY_);
    set(S_fig.hAx_mean, 'YLim', [0, 100]); % TW
    xylabel_(S_fig.hAx_mean, 'Time (sec)', sprintf('Common ref. (MAD, %s)', S_fig.vcRef_view));
    fThresh_ref = strcmpi(S_fig.vcRef_view, 'binned') && ~isempty(S_fig.blank_thresh) && S_fig.blank_thresh ~= 0;
    if fThresh_ref
        set(S_fig.hPlot_mean_thresh, 'XData', vrTime_sec([1,end]), 'YData', repmat(S_fig.blank_thresh, [1,2]));
    else
        clearPlots(S_fig.hPlot_mean_thresh);
    end

    %-----
    % Traces plot
    if isempty(get_(S_fig, 'hPlot_traces'))
        S_fig.hPlot_traces = plot(S_fig.hAx_traces, nan, nan, 'Color', [1 1 1]*.5);
        S_fig.hPlot_traces_spk = plot(S_fig.hAx_traces, nan, nan, 'm.-', 'LineWidth', 1.5);
        S_fig.hPlot_traces_spk1 = plot(S_fig.hAx_traces, nan, nan, 'ro');
        S_fig.hPlot_traces_thresh = plot(S_fig.hAx_traces, nan, nan, 'm:');
        S_fig.hPlot_traces_bad = plot(S_fig.hAx_traces, nan, nan, 'r');
    end
    if S_fig.fFilter
        vcFilter = S_fig.vcFilter;
        mrWav_ = bit2uV_(S_fig.mnWav_filt(viPlot,:), struct_add_(P, vcFilter));
    else
        mrWav_ = meanSubtract(single(S_fig.mnWav_clean(viPlot,:)) * P.uV_per_bit);
    end
    multiplot(S_fig.hPlot_traces, S_fig.maxAmp, vrTime_sec, mrWav_, 1:nSites);
    if ~isempty(S_fig.viSite_bad)
        multiplot(S_fig.hPlot_traces_bad, S_fig.maxAmp, vrTime_sec, mrWav_(:,S_fig.viSite_bad), S_fig.viSite_bad);
    else
        clearPlots(S_fig.hPlot_traces_bad);
    end
    if S_fig.fShow_spk
        vlSpk_ = S_fig.spikeTimes >= S_fig.nlim_bin(1) & S_fig.spikeTimes <= S_fig.nlim_bin(end);
        spikeTimes_ = single(S_fig.spikeTimes(vlSpk_)-S_fig.nlim_bin(1)+1);
        vrTime_spk_ = single(S_fig.spikeTimes(vlSpk_)) / P.sampleRateHz;
        spikeSites_ = single(S_fig.spikeSites(vlSpk_));
    else
        spikeTimes_ = [];
    end
    if isempty(spikeTimes_)
        clearPlots(S_fig.hPlot_traces_spk1);
        menu_label_('menu_preview_view_spike', 'Show [S]pikes');
    else
        multiplot(S_fig.hPlot_traces_spk1, S_fig.maxAmp, vrTime_spk_, mr2vr_sub2ind_(mrWav_, spikeTimes_, spikeSites_), spikeSites_, 1);
        menu_label_('menu_preview_view_spike', 'Hide [S]pikes');
    end
    vrThresh_site_uV = bit2uV_(-S_fig.siteThresholds(:), setfield(P, 'vcFilter', S_fig.vcFilter));
    vrThresh_site_uV(S_fig.viSite_bad) = nan;
    if S_fig.fThresh_spk && S_fig.fFilter
        multiplot(S_fig.hPlot_traces_thresh, S_fig.maxAmp, vrTime_sec([1,end,end])', repmat(vrThresh_site_uV, [1,3])');
        multiplot(S_fig.hPlot_traces_spk, S_fig.maxAmp, vrTime_sec, ...
        mr_set_(mrWav_, ~S_fig.mlWav_thresh(viPlot,:), nan)); %show spikes
        menu_label_('menu_preview_view_threshold', 'Hide spike [T]threshold');
    else
        clearPlots([S_fig.hPlot_traces_thresh, S_fig.hPlot_traces_spk]);
        menu_label_('menu_preview_view_threshold', 'Show spike [T]hreshold');
    end
    xylabel_(S_fig.hAx_traces, '', 'Site #');
    vcFilter_ = ifeq_(S_fig.fFilter, sprintf('Filter=%s', S_fig.vcFilter), 'Filter off');
    set(hFig, 'Name', sprintf('%s; %s; CommonRef=%s', P.paramFile, vcFilter_, S_fig.vcCommonRef));

    title_(S_fig.hAx_traces, sprintf('Scale: %0.1f uV', S_fig.maxAmp));
    menu_label_('menu_preview_view_filter', ifeq_(S_fig.fFilter, 'Show raw traces [F]', 'Show [F]iltered traces'));
    menu_label_('menu_preview_view_grid', ifeq_(S_fig.fGrid, 'Hide [G]rid', 'Show [G]rid'));
    if ~fKeepView
        set(S_fig.hAx_traces, 'YTick', 1:nSites, 'YLim', S_fig.siteLim + [-1,1], 'XLim', tlim_sec);
    end


    %-----
    % Site plot
    if isempty(get_(S_fig, 'hPlot_site'))
        S_fig.hPlot_site = barh(S_fig.hAx_sites, nan, nan, 1);
        S_fig.hPlot_site_bad = barh(S_fig.hAx_sites, nan, nan, 1, 'r');
        S_fig.hPlot_site_thresh = plot(S_fig.hAx_sites, nan, nan, 'r');
    end
    switch S_fig.vcSite_view
        case 'Site correlation', vrPlot_site = S_fig.vrCorr_max_site;
        case 'Spike threshold', vrPlot_site = single(S_fig.siteThresholds);
        case 'Event rate (Hz)', vrPlot_site = S_fig.vrEventRate_site;
        case 'Event SNR (median)', vrPlot_site = S_fig.vrEventSnr_site;
    end
    set(S_fig.hPlot_site, 'XData', 1:nSites, 'YData', vrPlot_site); %switch statement
    xylabel_(S_fig.hAx_sites, S_fig.vcSite_view, 'Site #');
    set(S_fig.hAx_sites, 'YLim', S_fig.siteLim + [-1,1]);
    if isempty(S_fig.thresh_corr_bad_site) || ~strcmpi(S_fig.vcSite_view, 'Site correlation')
        clearPlots(S_fig.hPlot_site_thresh);
    else
        set(S_fig.hPlot_site_thresh, 'XData', S_fig.thresh_corr_bad_site *[1,1], 'YData', [0, nSites+1]);
    end
    if ~isempty(S_fig.viSite_bad);
        vrPlot_site_bad = vrPlot_site;
        vrPlot_site_bad(~S_fig.vlSite_bad) = 0;
        set(S_fig.hPlot_site_bad, 'XData', 1:nSites, 'YData', vrPlot_site_bad); %switch statement
    else
        clearPlots(S_fig.hPlot_site_bad);
    end
    title_(S_fig.hAx_sites, sprintf('thresh_corr_bad_site=%0.4f', S_fig.thresh_corr_bad_site));


    %-----
    % PSD plot
    if isempty(get_(S_fig, 'hPlot_psd'))
        S_fig.hPlot_psd = plot(S_fig.hAx_psd, nan, nan, 'k');
        S_fig.hPlot_clean_psd = plot(S_fig.hAx_psd, nan, nan, 'g');
        S_fig.hPlot_psd_thresh = plot(S_fig.hAx_psd, nan, nan, 'r');
    end
    set(S_fig.hPlot_psd, 'XData', S_fig.vrFreq_psd, 'YData', S_fig.vrPower_psd);
    set(S_fig.hPlot_clean_psd, 'XData', S_fig.vrFreq_psd, 'YData', S_fig.vrPower_clean_psd);
    xylabel_(S_fig.hAx_psd, 'Frequency (Hz)', 'Power [dB]', 'TODO: before and after cleaning');
    set(S_fig.hAx_psd, 'XLim', [0, P.sampleRateHz/2]);
    title_(S_fig.hAx_psd, sprintf('fft_thresh=%s', num2str(S_fig.fft_thresh)));

    grid_([S_fig.hAx_traces, S_fig.hAx_mean, S_fig.hAx_sites, S_fig.hAx_psd], S_fig.fGrid);

    % Exit
    set(hFig, 'UserData', S_fig);
    figure_wait_(0, hFig);
    tryClose(hWait);
end % function
