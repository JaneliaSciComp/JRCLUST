% @TODO: run spike detection and show detected spikes, or below threshold values
%function [hFigPreview, hFigPreview.figData] = doPlotFigPreview(hCfg, fKeepView)
function hFigPreview = doPlotFigPreview(hFigPreview, hCfg, fKeepView)
    %DOPLOTFIGPREVIEW
    if nargin < 2
        fKeepView = 0;
    end
    hWait = jrclust.utils.qMsgBox('Plotting...', 0, 1);

    nSites = size(hFigPreview.figData.tracesFilt, 2);

    viPlot = hFigPreview.figData.windowBounds(1):hFigPreview.figData.windowBounds(2);
    XData = viPlot / hCfg.sampleRate;
    tlim_sec = (hFigPreview.figData.windowBounds + [-1 1]) / hCfg.sampleRate;

    %-----
    % Mean plot
    if ~hFigPreview.hasPlot('hPlotMean')
        hFigPreview.addPlot('hPlotMean', @plot, hFigPreview.hAxes('hAxMean'), nan, nan, 'k');
        hFigPreview.addPlot('hPlotMeanThresh', @plot, hFigPreview.hAxes('hAxMean'), nan, nan, 'r');
    end

    if strcmp(hFigPreview.figData.refView, 'original')
        YData = hFigPreview.figData.vrWav_filt_mean(viPlot); % TODO: vrWav_filt_mean is never set
    else % binned
        YData = hFigPreview.figData.channelMeansMAD(viPlot);
    end

    hFigPreview.plotApply('hPlotMean', @set, 'XData', XData, 'YData', YData);
    hFigPreview.axApply('hAxMean', @set, 'YLim', [0, 100]);
    hFigPreview.axApply('hAxMean', @xlabel, 'Time (sec)');
    hFigPreview.axApply('hAxMean', @ylabel, sprintf('Common ref. (MAD, %s)', hFigPreview.figData.refView));

    fThreshRef = strcmpi(hFigPreview.figData.refView, 'binned') && ~isempty(hFigPreview.figData.blankThresh) && hFigPreview.figData.blankThresh ~= 0;
    if fThreshRef
        hFigPreview.plotApply('hPlotMeanThresh', @set, 'XData', XData([1, end]), 'YData', repmat(hFigPreview.figData.blankThresh, [1, 2]));
    else
        hFigPreview.hidePlot('hPlotMeanThresh');
    end

    %-----
    % Traces plot
    if ~hFigPreview.hasPlot('hPlotTraces')
        hFigPreview.addPlot('hPlotTraces', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'Color', [1 1 1]*.5);
        hFigPreview.addPlot('hPlot_traces_spk', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'm.-', 'LineWidth', 1.5);
        hFigPreview.addPlot('hPlot_traces_spk1', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'ro');
        hFigPreview.addPlot('hPlotTracesThresh', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'm:');
        hFigPreview.addPlot('hPlotTracesBad', @plot, hFigPreview.hAxes('hAxTraces'), nan, nan, 'r');
    end

    if hFigPreview.figData.fFilter
        hCfg.setTemporaryParams('filterType', hFigPreview.figData.filterType);
        mrWav_ = jrclust.utils.bit2uV(hFigPreview.figData.tracesFilt(viPlot, :), hCfg);
        hCfg.resetTemporaryParams();
    else
        mrWav_ = jrclust.utils.meanSubtract(single(hFigPreview.figData.tracesClean(viPlot, :))*hCfg.bitScaling);
    end

    hFigPreview.multiplot('hPlotTraces', hFigPreview.figData.maxAmp, XData, mrWav_, 1:nSites);

    % plot bad sites in red
    if ~isempty(hFigPreview.figData.ignoreSites)
        hFigPreview.multiplot('hPlotTracesBad', hFigPreview.figData.maxAmp, XData, mrWav_(:, hFigPreview.figData.ignoreSites), hFigPreview.figData.ignoreSites);
    else
        hFigPreview.hidePlot('hPlotTracesBad');
    end

    if hFigPreview.figData.fShow_spk
        inBounds = hFigPreview.figData.spikeTimes >= hFigPreview.figData.windowBounds(1) & hFigPreview.figData.spikeTimes <= hFigPreview.figData.windowBounds(end);
        spikeTimes = single(hFigPreview.figData.spikeTimes(inBounds) - hFigPreview.figData.windowBounds(1)+1);
        spikeTimesSec = single(hFigPreview.figData.spikeTimes(inBounds)) / hCfg.sampleRate;
        spikeSites = single(hFigPreview.figData.spikeSites(inBounds));
    else
        spikeTimes = [];
    end

    if isempty(spikeTimes)
        hFigPreview.hidePlot('hPlot_traces_spk1');
        %menu_label_('menu_preview_view_spike', 'Show [S]pikes');
    else
        hFigPreview.multiplot('hPlot_traces_spk1', hFigPreview.figData.maxAmp, spikeTimesSec, jrclust.utils.rowColSelect(mrWav_, spikeTimes, spikeSites), spikeSites, 1);
        %menu_label_('menu_preview_view_spike', 'Hide [S]pikes');
    end

    hCfg.setTemporaryParams('filterType', hFigPreview.figData.filterType);
    siteThreshuV = jrclust.utils.bit2uV(-hFigPreview.figData.siteThresh(:), hCfg);
    hCfg.resetTemporaryParams();
    siteThreshuV(hFigPreview.figData.ignoreSites) = nan;

    if hFigPreview.figData.fShowThresh && hFigPreview.figData.fFilter
        hFigPreview.multiplot('hPlotTracesThresh', hFigPreview.figData.maxAmp, XData([1,end,end])', repmat(siteThreshuV, [1, 3])');
        hFigPreview.multiplot('hPlot_traces_spk', hFigPreview.figData.maxAmp, XData, ...
            mr_set_(mrWav_, ~hFigPreview.figData.isThreshCrossing(viPlot, :), nan)); % show spikes
        %menu_label_('menu_preview_view_threshold', 'Hide spike [T]threshold');
    else
        hFigPreview.hidePlot('hPlotTracesThresh');
        hFigPreview.hidePlot('hPlot_traces_spk');
        %menu_label_('menu_preview_view_threshold', 'Show spike [T]hreshold');
    end

    hFigPreview.axApply('hAxTraces', @ylabel, 'Site #');
    vcFilter_ = jrclust.utils.ifEq(hFigPreview.figData.fFilter, sprintf('Filter=%s', hFigPreview.figData.filterType), 'Filter off');
    set(hFigPreview, 'Name', sprintf('%s; %s; CommonRef=%s', hCfg.vcFile_prm, vcFilter_, hFigPreview.figData.vcCommonRef));

    title_(hFigPreview.figData.hAxTraces, sprintf('Scale: %0.1f uV', hFigPreview.figData.maxAmp));
    menu_label_('menu_preview_view_filter', jrclust.utils.ifEq(hFigPreview.figData.fFilter, 'Show raw traces [F]', 'Show [F]iltered traces'));
    menu_label_('menu_preview_view_grid', jrclust.utils.ifEq(hFigPreview.figData.fGrid, 'Hide [G]rid', 'Show [G]rid'));
    if ~fKeepView
        set(hFigPreview.figData.hAxTraces, 'YTick', 1:nSites, 'YLim', hFigPreview.figData.siteLim + [-1,1], 'XLim', tlim_sec);
    end


    %-----
    % Site plot
    if isempty(get_(hFigPreview.figData, 'hPlot_site'))
        hFigPreview.figData.hPlot_site = barh(hFigPreview.figData.hAx_sites, nan, nan, 1);
        hFigPreview.figData.hPlot_site_bad = barh(hFigPreview.figData.hAx_sites, nan, nan, 1, 'r');
        hFigPreview.figData.hPlot_site_thresh = plot(hFigPreview.figData.hAx_sites, nan, nan, 'r');
    end
    switch hFigPreview.figData.siteView
        case 'Site correlation', vrPlot_site = hFigPreview.figData.maxCorrSite;
        case 'Spike threshold', vrPlot_site = single(hFigPreview.figData.siteThresh);
        case 'Event rate (Hz)', vrPlot_site = hFigPreview.figData.siteEventRate;
        case 'Event SNR (median)', vrPlot_site = hFigPreview.figData.siteEventSNR;
    end
    set(hFigPreview.figData.hPlot_site, 'XData', 1:nSites, 'YData', vrPlot_site); %switch statement
    xylabel_(hFigPreview.figData.hAx_sites, hFigPreview.figData.siteView, 'Site #');
    set(hFigPreview.figData.hAx_sites, 'YLim', hFigPreview.figData.siteLim + [-1,1]);
    if isempty(hFigPreview.figData.thresh_corr_bad_site) || ~strcmpi(hFigPreview.figData.siteView, 'Site correlation')
        hide_plot_(hFigPreview.figData.hPlot_site_thresh);
    else
        set(hFigPreview.figData.hPlot_site_thresh, 'XData', hFigPreview.figData.thresh_corr_bad_site *[1,1], 'YData', [0, nSites+1]);
    end
    if ~isempty(hFigPreview.figData.ignoreSites);
        vrPlot_site_bad = vrPlot_site;
        vrPlot_site_bad(~hFigPreview.figData.vlSite_bad) = 0;
        set(hFigPreview.figData.hPlot_site_bad, 'XData', 1:nSites, 'YData', vrPlot_site_bad); %switch statement
    else
        hide_plot_(hFigPreview.figData.hPlot_site_bad);
    end
    title_(hFigPreview.figData.hAx_sites, sprintf('thresh_corr_bad_site=%0.4f', hFigPreview.figData.thresh_corr_bad_site));


    %-----
    % PSD plot
    if isempty(get_(hFigPreview.figData, 'hPlot_psd'))
        hFigPreview.figData.hPlot_psd = plot(hFigPreview.figData.hAxPSD, nan, nan, 'k');
        hFigPreview.figData.hPlot_clean_psd = plot(hFigPreview.figData.hAxPSD, nan, nan, 'g');
        hFigPreview.figData.hPlot_psd_thresh = plot(hFigPreview.figData.hAxPSD, nan, nan, 'r');
    end
    set(hFigPreview.figData.hPlot_psd, 'XData', hFigPreview.figData.vrFreq_psd, 'YData', hFigPreview.figData.vrPower_psd);
    set(hFigPreview.figData.hPlot_clean_psd, 'XData', hFigPreview.figData.vrFreq_psd, 'YData', hFigPreview.figData.vrPower_clean_psd);
    xylabel_(hFigPreview.figData.hAxPSD, 'Frequency (Hz)', 'Power [dB]', 'TODO: before and after cleaning');
    set(hFigPreview.figData.hAxPSD, 'XLim', [0, hCfg.sampleRate/2]);
    title_(hFigPreview.figData.hAxPSD, sprintf('fft_thresh=%s', num2str(hFigPreview.figData.fft_thresh)));

    grid_([hFigPreview.figData.hAxTraces, hFigPreview.figData.hAxMean, hFigPreview.figData.hAx_sites, hFigPreview.figData.hAxPSD], hFigPreview.figData.fGrid);

    % Exit
    set(hFigPreview, 'UserData', hFigPreview.figData);
    figure_wait_(0, hFigPreview);
    close_(hWait);
end %func

%% LOCAL FUNCTIONS
function mr = mr_set_(mr, ml, val)
    mr(ml)=val;
end %func
