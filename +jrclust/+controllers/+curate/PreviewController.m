classdef PreviewController < handle
    %PREVIEWCONTROLLER
    
    properties (SetAccess=private)
        fileBounds;
        hCfg;
        hFigPreview;
        tracesRaw;
    end

    properties (Hidden, SetAccess=private)
        hAxPSD;
        hAxSites;
        hAxTraces;
        nLoads;
        nLoadsPerFile;
        nSamplesPerLoad;
    end

    %% LIFECYCLE
    methods
        function obj = PreviewController(hCfg)
            %PREVIEWCONTROLLER Construct an instance of this class
            obj.hCfg = hCfg;
            obj.fileBounds = containers.Map();
        end
    end

    %% KEYBOARD/MOUSE FUNCTIONS
    methods (Hidden)
        function keyPressFigPreview(obj, hObject, hEvent)
            S_fig = obj.hFigPreview.figData;
            factor = 1 + 3 * keyMod(hEvent, 'shift');
            nSites = numel(obj.hCfg.viSite2Chan);

            switch lower(hEvent.Key)
                case {'uparrow', 'downarrow'}
                    S_fig.maxAmp = change_amp_(hEvent, S_fig.maxAmp, ...
                        S_fig.hPlotTraces, S_fig.hPlotTracesBad, S_fig.hPlotTracesThresh, ...
                        S_fig.hPlot_traces_spk, S_fig.hPlot_traces_spk1);
                    title_(obj.hAxTraces, sprintf('Scale: %0.1f uV', S_fig.maxAmp));
                    set(hObject, 'UserData', S_fig);

                case {'leftarrow', 'rightarrow', 'home', 'end'} %navigation
                    switch lower(hEvent.Key)
                        case 'leftarrow'
                            windowBounds = S_fig.windowBounds - (S_fig.windowWidth) * factor; %no overlap
                            if windowBounds(1)<1
                                jrclust.utils.qMsgBox('Beginning of file', 1);
                                windowBounds = [1, S_fig.windowWidth];
                            end
                        case 'rightarrow'
                            windowBounds = S_fig.windowBounds + (S_fig.windowWidth + 1) * factor; %no overlap
                            if windowBounds(2) > S_fig.nSamplesTotal
                                jrclust.utils.qMsgBox('End of file', 1);
                                windowBounds = [-S_fig.windowWidth+1, 0] + S_fig.nSamplesTotal;
                            end
                        case 'home', windowBounds = [1, S_fig.windowWidth]; %begging of file
                        case 'end', windowBounds = [-S_fig.windowWidth+1, 0] + S_fig.nSamplesTotal; %end of file
                    end %switch
                    S_fig.windowBounds = windowBounds;
                    set(hObject, 'UserData', S_fig);
                    doPlotFigPreview(hCfg);

                case 'f' %apply filter
                    set(hObject, 'UserData', setfield(S_fig, 'fFilter', ~S_fig.fFilter));
                    doPlotFigPreview(hCfg, 1);

                case 't' %toggle spike threshold
                    set(hObject, 'UserData', setfield(S_fig, 'fShowThresh', ~S_fig.fShowThresh));
                    doPlotFigPreview(hCfg, 1);

                case 's' %show/hide spikes
                    set(hObject, 'UserData', setfield(S_fig, 'fShow_spk', ~S_fig.fShow_spk));
                    doPlotFigPreview(hCfg, 1);

                case 'g' %grid toggle on/off
                    S_fig.fGrid = ~S_fig.fGrid;
                    grid_([obj.hAxTraces, S_fig.hAxMean, S_fig.hAxPSD, S_fig.hAx_sites], S_fig.fGrid);
                    set(hObject, 'UserData', S_fig);
                    menu_label_('menu_preview_view_grid', jrclust.utils.ifEq(S_fig.fGrid, 'Hide [G]rid', 'Show [G]rid'));

                case 'h' % help
                    jrclust.utils.qMsgBox(S_fig.helpText, 1);

                case 'r' % reset view
                    axis_(obj.hAxTraces, [S_fig.windowBounds / hCfg.sRateHz, S_fig.siteLim+[-1,1]]);

                case 'e' % export to workspace
                    Fig_preview_export_(hObject);

            end %switch
        end %func
    end

    %% UTILITY METHODS
    methods (Hidden)
        function addMenu(obj)
            %ADDMENU
            obj.hFigPreview.figApply(@set, 'MenuBar', 'None');
            fileMenu = obj.hFigPreview.figApply(@uimenu, 'Label', 'File');
            uimenu(fileMenu, 'Label', sprintf('Save to %s', obj.hCfg.configFile), 'Callback', @(hO, hE) obj.savePrm());
            %%%%%%%%%%%%%%%%
            uimenu(fileMenu, 'Label', '[E]xport to workspace', 'Callback', @(hO, hE) Fig_preview_export_(hFig));
            uimenu(fileMenu, 'Label', 'Save spike detection threshold', 'Callback', @(hO, hE) Fig_preview_save_threshold_(hFig));

            % Edit menu
            mh_edit = uimenu(hFig, 'Label','Edit');

            uimenu(mh_edit, 'Label', 'Bad site threshold', 'Callback', @(hO, hE) Fig_preview_site_thresh_(hFig));
            uimenu(mh_edit, 'Label', 'Spike detection threshold', 'Callback', @(hO, hE) Fig_preview_spk_thresh_(hFig));

            mh_edit_filter = uimenu(mh_edit, 'Label', 'Filter mode');
            uimenu_options_(mh_edit_filter, {'ndiff', 'bandpass', 'sgdiff', 'fir1', 'user'}, @Fig_preview_filter_, hFig);
            menu_checkbox_(mh_edit_filter, get_filter_(obj.hCfg));

            mh_edit_ref = uimenu(mh_edit, 'Label', 'Reference mode');
            uimenu_options_(mh_edit_ref, {'none', 'mean', 'median'}, @Fig_preview_ref_, hFig); % @TODO: local mean
            menu_checkbox_(mh_edit_ref, obj.hCfg.vcCommonRef);

            uimenu(mh_edit, 'Label', 'Common reference threshold', 'Callback', @(hO, hE)Fig_preview_ref_thresh_(hFig));

            uimenu(mh_edit, 'Label', 'FFT cleanup threshold', 'Callback', @(hO, hE)Fig_preview_fft_thresh_(hFig));


            % View menu
            mh_view = uimenu(hFig, 'Label','View');

            mh_view_trange = uimenu(mh_view, 'Label', 'Display time range (s)');
            uimenu_options_(mh_view_trange, {'0.05', '0.1', '0.2', '0.5', '1', '2', '5', 'Custom'}, @Fig_preview_trange_, hFig);
            menu_checkbox_(mh_view_trange, '0.1');
            uimenu(mh_view, 'Label', 'Display site range', 'Callback', @(hO, hE)Fig_preview_site_range_(hFig));

            uimenu(mh_view, 'Label', 'Show raw traces [F]', 'Callback', @(hO, hE)keyPressFcn_cell_(hFig, 'f'), 'Tag', 'menu_preview_view_filter');
            uimenu(mh_view, 'Label', 'Show spike [T]hreshold', 'Callback', @(hO, hE)keyPressFcn_cell_(hFig, 't'), 'Tag', 'menu_preview_view_threshold');
            uimenu(mh_view, 'Label', 'Hide [S]pikes', 'Callback', @(hO, hE)keyPressFcn_cell_(hFig, 's'), 'Tag', 'menu_preview_view_spike');
            uimenu(mh_view, 'Label', 'Hide [G]rid', 'Callback', @(hO, hE)keyPressFcn_cell_(hFig, 'g'), 'Tag', 'menu_preview_view_grid');

            mh_view_site = uimenu(mh_view, 'Label', 'Site view');
            uimenu_options_(mh_view_site, {'Site correlation', 'Spike threshold', 'Event rate (Hz)', 'Event SNR (median)'}, @Fig_preview_site_plot_, hFig);
            menu_checkbox_(mh_view_site, 'Site correlation');

            mh_view_ref = uimenu(mh_view, 'Label', 'Reference view');
            uimenu_options_(mh_view_ref, {'original', 'binned'}, @Fig_preview_view_ref_, hFig);
            menu_checkbox_(mh_view_ref, 'original');

            mh_view_freq = uimenu(mh_view, 'Label', 'Frequency scale');
            uimenu_options_(mh_view_freq, {'Linear', 'Log'}, @Fig_preview_psd_plot_, hFig);
            menu_checkbox_(mh_view_freq, 'Linear');

            % mh_view_psd = uimenu(mh_view, 'Label', 'PSD view');
            % uimenu_options_(mh_view_psd, {'original', 'detrended'}, @Fig_preview_view_psd_, hFig);
            % menu_checkbox_(mh_view_psd, 'original');

            % 'Power', 'Detrended',
        end %func

        function loadPreview(obj)
            %LOADPREVIEW
            nLoadsMax = obj.hCfg.nLoads_max_preview;
%             nSecsLoad = obj.hCfg.sec_per_load_preview;
            nSecsLoad = 1;

            rawRecordings = jrclust.utils.subsample(obj.hCfg.rawRecordings, nLoadsMax);

            % load files
            nFiles = numel(rawRecordings);
            obj.nLoadsPerFile = floor(nLoadsMax / nFiles);
            obj.nSamplesPerLoad = ceil(nSecsLoad*obj.hCfg.sampleRate);

            obj.hCfg.useGPU = false;

            tracesRaw_ = cell(nFiles, 1);

            for iFile = 1:nFiles
                hRec = jrclust.models.recording.Recording(rawRecordings{iFile}, obj.hCfg);
                if hRec.nSamples <= obj.nSamplesPerLoad
                    nLoadsFile = 1;
                    nSamplesLoad = hRec.nSamples;
                else
                    nLoadsFile = min(obj.nLoadsPerFile, floor(hRec.nSamples / obj.nSamplesPerLoad));
                    nSamplesLoad = obj.nSamplesPerLoad;
                end

                multiBounds = sample_skip_([1, nSamplesLoad], hRec.nSamples, nLoadsFile);
                obj.fileBounds(hRec.binpath) = multiBounds;

                fileTraces_ = cell(nLoadsFile, 1);

                for iLoad = 1:nLoadsFile
                    iBounds = multiBounds{iLoad};
                    fileTraces_{iLoad} = hRec.readROI(obj.hCfg.siteMap, iBounds(1):iBounds(2))';
                end

                tracesRaw_{iFile} = cat(1, fileTraces_{:});
            end

            obj.tracesRaw = cat(1, tracesRaw_{:});
        end

        function savePrm(obj)
            %SAVEPRM
            % Update to a parameter file, show a preview window
            % List of parameters to update. Set threshold: save to a known location and load
            if isempty(obj.hFigPreview.figData)
                return;
            end

            % Select variables to export
            newPrms = struct('fft_thresh', obj.hFigPreview.figData.fft_thresh, ...
                             'viSiteZero', find(obj.hFigPreview.figData.vlSite_bad), ...
                             'qqFactor', obj.hFigPreview.figData.qqFactor, ...
                             'filterType', obj.hFigPreview.figData.filterType, ...
                             'vcCommonRef', obj.hFigPreview.figData.vcCommonRef, ...
                             'blankThresh', obj.hFigPreview.figData.blankThresh, ...
                             'blank_period_ms', obj.hFigPreview.figData.blank_period_ms);

            % preview and edit variables in the edit box
            strNewPrms = jrclust.utils.struct2str(newPrms);
            dlgAns = inputdlg(obj.hCfg.configFile, 'Update confirmation', 16, {strNewPrms}, struct('Resize', 'on'));

            if isempty(dlgAns)
                return;
            end

            newPrms = jrclust.utils.str2struct(dlgAns{1});
            if isempty(newPrms)
                return;
            end

            obj.hCfg.fftThreshMAD = newPrms.fft_thresh;
            obj.hCfg.ignoreSites = newPrms.viSiteZero;
            obj.hCfg.qqFactor = newPrms.qqFactor;
            obj.hCfg.filterType = newPrms.filterType;
            obj.hCfg.carMode = newPrms.vcCommonRef;
            obj.hCfg.blankThresh = newPrms.blankThresh;
            obj.hCfg.blank_period_ms = newPrms.blank_period_ms;

            obj.hCfg.flush();
            obj.hCfg.edit();
        end

        function updateFigPreview(obj)
        end
    end

    %% USER METHODS
    methods
        function preview(obj)
            obj.loadPreview();

            nSites = size(obj.tracesRaw, 2);

            % process signal, how about common mean?

            siteCorr = corr(single(obj.tracesRaw));
            siteCorr(logical(eye(size(siteCorr)))) = 0;

            % Create a Figure
            gap = .05;
            obj.hFigPreview = jrclust.views.Figure('FigPreview', [0 0 .5 1], obj.hCfg.configFile, true, true); %plot a summary pannel
            obj.hFigPreview.addAxes('hAxMean', 'Position', [gap, gap, 3/4-gap, 1/4-gap], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxTraces', 'Position', [gap, 1/4+gap, 3/4-gap, 3/4-gap*2], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxSites', 'Position', [3/4+gap, gap, 1/4-gap*1.5, 2/3-gap*2], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxPSD', 'Position', [3/4+gap, 2/3+gap, 1/4-gap*1.5, 1/3-gap*2], 'NextPlot', 'add');
            linkaxes([obj.hFigPreview.hAxes('hAxMean'), obj.hFigPreview.hAxes('hAxTraces')], 'x');

            % Callback functions
            %obj.addMenu();
            obj.hFigPreview.figApply(@set, 'KeyPressFcn', @obj.keyPressFigPreview, 'BusyAction', 'cancel');
            %mouse_figure(obj.hFigPreview, hAxTraces);
            obj.hFigPreview.setMouseable([], 'hAxTraces');

            % Build figData
            figData = struct('filterType', obj.hCfg.filterType, ...
                             'carMode', obj.hCfg.carMode, ...
                             'siteCorrThresh', obj.hCfg.siteCorrThresh, ...
                             'fftThreshMAD', obj.hCfg.fftThreshMAD, ...
                             'qqFactor', obj.hCfg.qqFactor, ...
                             'blankThresh', obj.hCfg.blankThresh, ...
                             'blank_period_ms', obj.hCfg.blank_period_ms, ...
                             'ignoreSites', obj.hCfg.ignoreSites, ...
                             'nSamplesTotal', size(obj.tracesRaw, 1), ...
                             'maxAmp', obj.hCfg.maxAmp, ...
                             'nLoads', obj.nLoads, ...
                             'tracesRaw', obj.tracesRaw, ...
                             'maxCorrSite', max(siteCorr), ...
                             'siteLim', [1 nSites], ...
                             'fFilter', true, ...
                             'fGrid', true, ...
                             'fShowThresh', false, ...
                             'fShow_spk', true, ...
                             'siteView', 'Site correlation', ...
                             'refView', 'binned', ...
                             'psdView', 'original');

            if isprop(obj.hCfg, 'preview_window')
                figData.windowWidth = round(obj.hCfg.preview_window * obj.hCfg.sampleRate); % TW
            else
                figData.windowWidth = obj.nSamplesPerLoad;
            end

            figData.windowBounds = [1, figData.windowWidth];
            figData.helpText = {'Left/Right: change time (Shift: x4)', ...
                                '[Home/End]: go to beginning/end of file', ...
                                '---------', ...
                                'Up/Down: change scale (Shift: x4)', ...
                                'Right/Left: change time (Shift: x4)', ...
                                'Zoom: Mouse wheel', ...
                                '[x/y/ESC]: zoom direction', ...
                                'Pan: hold down the wheel and drag', ...
                                '---------', ...
                                '[F]ilter toggle', ...
                                '[G]rid toggle'};

            obj.hFigPreview.figData = figData;

            drawnow;

            doUpdateFigPreview(obj.hFigPreview, figData, false, obj.hCfg);
        end
    end
end

