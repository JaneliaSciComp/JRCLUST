%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function S_fig = preview_(P, fDebug_ui_)
    % Data summary figure, interactive

    global fDebug_ui
    if nargin<2, fDebug_ui_ = 0; end
    fDebug_ui = fDebug_ui_;
    setUserData(fDebug_ui);
    if ischar(P), P = loadParams(P); end
    [mnWav_raw, S_preview] = load_preview_(P);
    setUserData(P);
    nSites = size(mnWav_raw,2);

    % process signal, how about common mean?

    % Bad channel metrics, do it once
    mrCorr_site = corr(single(mnWav_raw));
    mrCorr_site(logical(eye(size(mrCorr_site)))) = 0;
    vrCorr_max_site = max(mrCorr_site);

    % Create a Figure
    gap = .05;
    hFig = createFigure('Fig_preview', [0 0 .5 1], P.paramFile, 1, 1); %plot a summary pannel
    hAx_mean = axes('Parent', hFig, 'Position',      [gap        gap         3/4-gap         1/4-gap], 'NextPlot', 'add');
    hAx_traces = axes('Parent', hFig, 'Position',    [gap        1/4+gap     3/4-gap         3/4-gap*2], 'NextPlot', 'add');
    hAx_sites = axes('Parent', hFig, 'Position',     [3/4+gap,   0+gap       1/4-gap*1.5     2/3-gap*2], 'NextPlot', 'add');
    hAx_psd = axes('Parent', hFig, 'Position',       [3/4+gap,   2/3+gap     1/4-gap*1.5     1/3-gap*2], 'NextPlot', 'add');
    linkaxes([hAx_mean, hAx_traces], 'x');

    % Callback functions
    Fig_preview_menu_(hFig);
    set(hFig, 'KeyPressFcn', @keyPressFcn_Fig_preview_, 'BusyAction', 'cancel');
    mouse_figure(hFig, hAx_traces);

    % Build S_fig
    [nLoads, nSamples_bin, maxAmp] = deal(S_preview.nLoads, size(mnWav_raw,1), P.maxAmp);
    % nLoad_bin = S_preview.nSamples_per_load;
    nLoad_bin = round(P.preview_window * P.sampleRateHz); % TW
    nlim_bin = [1, nLoad_bin];
    siteLim = [1, nSites];
    [vcFilter, vcCommonRef, thresh_corr_bad_site, fft_thresh, qqFactor, blank_thresh, blank_period_ms, viSiteZero] = ...
    get_(P, 'vcFilter', 'vcCommonRef', 'thresh_corr_bad_site', 'fft_thresh', 'qqFactor', 'blank_thresh', 'blank_period_ms', 'viSiteZero');
    [fGrid, fFilter, fThresh_spk, fShow_spk] = deal(1, 1, 0, 1);
    [vcSite_view, vcRef_view, vcPsd_view] = deal('Site correlation', 'binned', 'original');
    csHelp = { ...
    'Left/Right: change time (Shift: x4)', ...
    '[Home/End]: go to beginning/end of file', ...
    '---------', ...
    'Up/Down: change scale (Shift: x4)', ...
    'Right/Left: change time (Shift: x4)', ...
    'Zoom: Mouse wheel', ...
    '[x/y/ESC]: zoom direction', ...
    'Pan: hold down the wheel and drag', ...
    '---------', ...
    '[F]ilter toggle', ...
    'Gri[D] toggle', ...
    };
    S_fig = makeStruct_(...
    vcFilter, vcCommonRef, thresh_corr_bad_site, fft_thresh, qqFactor, blank_thresh, blank_period_ms, viSiteZero, ...
    nlim_bin, nLoad_bin, nLoads, nSamples_bin, maxAmp, ...
    mnWav_raw, vrCorr_max_site, S_preview, csHelp, ...
    hAx_mean, hAx_traces, hAx_sites, hAx_psd, ...
    fFilter, fGrid, fThresh_spk, siteLim, fShow_spk, ...
    vcSite_view, vcRef_view, vcPsd_view);

    % Exit
    set(hFig, 'UserData', S_fig);
    drawnow;
    S_fig = Fig_preview_update_(hFig, S_fig, 0);
end %func
