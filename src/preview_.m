%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function S_fig = preview_(P, fDebug_ui_)
    % Data summary figure, interactive

    global fDebug_ui
    if nargin < 2
        fDebug_ui_ = 0;
    end
    fDebug_ui = fDebug_ui_;
    set0_(fDebug_ui);

    if ischar(P)
        P = loadParam_(P);
    end

    [mnWav_raw, S_preview] = load_preview_(P);
    set0_(P);
    nSites = size(mnWav_raw,2);

    % process signal, how about common mean?

    % Bad channel metrics, do it once
    mrCorr_site = corr(single(mnWav_raw));
    mrCorr_site(logical(eye(size(mrCorr_site)))) = 0;
    vrCorr_max_site = max(mrCorr_site);

    % Create a Figure
    gap = .05;
    hFig = create_figure_('Fig_preview', [0 0 .5 1], P.vcFile_prm, 1, 1); %plot a summary pannel
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
    % nLoad_bin = round(0.1 * P.sRateHz);
    if isfield(P, 'preview_window')
        nLoad_bin = round(P.preview_window * P.sampleRateHz); % TW
    else
        nLoad_bin = S_preview.nSamples_per_load;
    end

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

%% local functions
%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function [mnWav_raw, S_preview] = load_preview_(P)
    % Load the subsampled dataset
    % Useful for inspecting threshold and so on. filter and
    % S_preview: which file and where it came from
    if ischar(P)
        P = loadParam_(P);
    end

    nLoads_max_preview = get_set_(P, 'nLoads_max_preview', 30);
    sec_per_load_preview = get_set_(P, 'sec_per_load_preview', 1);

    % determine files to load
    if isempty(P.csFile_merge)
        csFile_bin = {P.vcFile};
    else
        csFile_bin = filter_files_(P.csFile_merge);
    end
    csFile_bin = subsample_(csFile_bin, nLoads_max_preview);

    % load files
    nLoads_per_file = floor(nLoads_max_preview / numel(csFile_bin));
    nSamples_per_load = round(sec_per_load_preview * P.sRateHz);

    % file loading loop
    [mnWav_raw, cviLim_load, csFile_load] = deal({});
    % [mnWav_raw, mnWav_filt] = deal({});
    P.fGpu = 0;
    for iFile = 1:numel(csFile_bin)
        try
            vcFile_bin_ = csFile_bin{iFile};
            [fid_bin_, nBytes_bin_, P.header_offset] = fopen_(vcFile_bin_, 'r');
            set0_(P);
            if isempty(fid_bin_)
                fprintf(2, '.bin file does not exist: %s\n', vcFile_bin_);
                continue;
            end
            fprintf('File %d/%d: %s\n\t', iFile, numel(csFile_bin), vcFile_bin_);
            nSamples_bin_ = floor(nBytes_bin_ / bytesPerSample_(P.vcDataType) / P.nChans);
            if nSamples_bin_ < nSamples_per_load % load the whole thing
                nLoads_per_file_ = 1;
                nSamples_per_load_ = nSamples_bin_;
            else
                nLoads_per_file_ = min(nLoads_per_file, floor(nSamples_bin_ / nSamples_per_load));
                nSamples_per_load_ = nSamples_per_load;
            end
            [cvi_lim_bin, viRange_bin] = sample_skip_([1, nSamples_per_load_], nSamples_bin_, nLoads_per_file_);
            for iLoad_ = 1:nLoads_per_file_
                fprintf('.');
                ilim_bin_ = cvi_lim_bin{iLoad_};

                % fseek_(fid_bin_, ilim_bin_(1), P);
                dshape_ = [P.nChans, diff(ilim_bin_) + 1];
                offset_ = max(0, (iSample_bin-1) * P.nChans * bytesPerSample_(P.vcDataType) + get_set_(P, 'header_offset', 0));

                % mnWav_raw{end+1} = load_file_(fid_bin_, diff(ilim_bin_) + 1, P);
                mnWav_raw{end+1} = jrclust.utils.readRecording(fid_bin_, P.vcDataType, dshape_, offset_, P);
                cviLim_load{end+1} = ilim_bin_;
                csFile_load{end+1} = vcFile_bin_;
            end
            fprintf('\n');
        catch
            fprintf(2, 'Loading error: %s\n', vcFile_bin_);
        end
        fclose_(fid_bin_, 0);
    end
    nLoads = numel(mnWav_raw);
    mnWav_raw = cell2mat(mnWav_raw');
    % if nargout>=2, mnWav_raw = cell2mat(mnWav_raw'); end
    if nargout>=2
        S_preview = makeStruct_(nLoads_per_file, nLoads_max_preview, ...
        sec_per_load_preview, nSamples_per_load, nLoads, csFile_load, cviLim_load);
    end
end %func

