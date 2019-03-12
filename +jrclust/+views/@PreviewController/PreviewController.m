classdef PreviewController < jrclust.interfaces.FigureController
    %PREVIEWCONTROLLER
    properties (SetAccess=private)
        figData;
        fileBounds;
        hCfg;
        hFigPreview;
    end

    properties (Dependent)
        blankPeriod;
        blankThresh;
        CARMode;
        channelMeansMAD;
        fFilter;
        fftThresh;
        fGrid;
        filterType;
        fShowSpikes;
        fShowThresh;
        helpText;
        ignoreMe;
        ignoreSites;
        isThreshCrossing;
        keepMe;
        maxAmp;
        maxCorrSite;
        nLoads;
        nSamplesTotal;
        psdFreq;
        psdPower;
        psdPowerClean;
        psdView;
        qqFactor;
        refView;
        siteCorrThresh;
        siteEventRate;
        siteEventSNR
        siteSNR;
        siteThresh;
        siteLim;
        siteView;
        spikeAmps;
        spikeSites;
        spikeTimes;
        tracesCAR;
        tracesClean;
        tracesFilt;
        tracesRaw;
        windowBounds;
        windowWidth;
    end

    properties (Hidden, SetAccess=private)
        nLoadsPerFile;
        nSamplesPerLoad;
    end

    %% LIFECYCLE
    methods
        function obj = PreviewController(hCfg)
            %PREVIEWCONTROLLER Construct an instance of this class
            obj.hCfg = hCfg;
            obj.fileBounds = containers.Map();
            obj.figData = struct();
        end
    end

    %% KEYBOARD/MOUSE FUNCTIONS
    methods (Hidden)
        function keyPressFigPreview(obj, ~, hEvent)
            factor = 4^double(jrclust.utils.keyMod(hEvent, 'shift')); % 1 or 4

            switch lower(hEvent.Key)
                case 'uparrow'
                    obj.maxAmp = obj.maxAmp*sqrt(2)^-factor;
                    obj.hFigPreview.rescalePlot('hPlotTraces', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlotTracesBad', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlotTracesThresh', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlot_traces_spk', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlot_traces_spk1', obj.maxAmp);

                    obj.hFigPreview.axApply('hAxTraces', @title, sprintf('Scale: %0.1f uV', obj.maxAmp));

                case 'downarrow'
                    obj.maxAmp = obj.maxAmp*sqrt(2)^factor;
                    obj.hFigPreview.rescalePlot('hPlotTraces', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlotTracesBad', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlotTracesThresh', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlot_traces_spk', obj.maxAmp);
                    obj.hFigPreview.rescalePlot('hPlot_traces_spk1', obj.maxAmp);

                    obj.hFigPreview.axApply('hAxTraces', @title, sprintf('Scale: %0.1f uV', obj.maxAmp));

                case {'leftarrow', 'rightarrow', 'home', 'end'}
                    switch lower(hEvent.Key)
                        case 'leftarrow'
                            if obj.windowBounds(1) == 1
                                jrclust.utils.qMsgBox('Beginning of file', 1);
                                return;
                            end

                            newBounds = obj.windowBounds - obj.windowWidth*factor;
                            if newBounds(1) < 1
                                newBounds = [1, obj.windowWidth];
                            end

                        case 'rightarrow'
                            if obj.windowBounds(end) == obj.nSamplesTotal
                                jrclust.utils.qMsgBox('End of file', 1);
                                return;
                            end

                            newBounds = obj.windowBounds + (obj.windowWidth + 1)*factor;
                            if newBounds(2) > obj.nSamplesTotal
                                newBounds = [1 - obj.windowWidth, 0] + obj.nSamplesTotal;
                            end

                        case 'home'
                            newBounds = [1, obj.windowWidth]; % beginning of file

                        case 'end'
                            newBounds = [1 - obj.windowWidth, 0] + obj.nSamplesTotal; % end of file
                    end % switch

                    if ~all(newBounds == obj.windowBounds) % don't update if there's no change
                        obj.windowBounds = newBounds;
                        doPlotFigPreview(obj.hFigPreview, obj.figData, 0, obj.hCfg);
                    end

                case 'e' % export to workspace
                    obj.exportTraces();

                case 'f' % apply filter
                    obj.fFilter = ~obj.fFilter;
                    doPlotFigPreview(obj.hFigPreview, obj.figData, 1, obj.hCfg);

                case 'g' % grid toggle on/off
                    obj.fGrid = ~obj.fGrid;
                    obj.hFigPreview.axApply('hAxMean', @grid, jrclust.utils.ifEq(obj.fGrid, 'on', 'off'));
                    obj.hFigPreview.axApply('hAxTraces', @grid, jrclust.utils.ifEq(obj.fGrid, 'on', 'off'));
                    obj.hFigPreview.axApply('hAxSites', @grid, jrclust.utils.ifEq(obj.fGrid, 'on', 'off'));
                    obj.hFigPreview.axApply('hAxPSD', @grid, jrclust.utils.ifEq(obj.fGrid, 'on', 'off'));

                case 'h' % help
                    jrclust.utils.qMsgBox(obj.helpText, 1);

                case 'r' % reset view
                    obj.hFigPreview.axApply('hAxTraces', @axis, [obj.windowBounds/obj.hCfg.sampleRate, obj.siteLim+[-1, 1]]);

                case 's' % show/hide spikes
                    obj.fShowSpikes = ~obj.fShowSpikes;
                    doPlotFigPreview(obj.hFigPreview, obj.figData, 1, obj.hCfg);

                case 't' % toggle spike threshold
                    obj.fShowThresh = ~obj.fShowThresh;
                    doPlotFigPreview(obj.hFigPreview, obj.figData, 1, obj.hCfg);

            end %switch
        end %func
    end

    %% UTILITY METHODS
    methods (Hidden)
        function addMenu(obj)
            %ADDMENU
            obj.hFigPreview.figApply(@set, 'MenuBar', 'None');

            % File menu
            fileMenu = obj.hFigPreview.figApply(@uimenu, 'Label', 'File');
            uimenu(fileMenu, 'Label', sprintf('Save to %s', obj.hCfg.configFile), 'Callback', @(hO, hE) obj.saveConfig());
            uimenu(fileMenu, 'Label', '[E]xport to workspace', 'Callback', @(hO, hE) obj.exportTraces());
            uimenu(fileMenu, 'Label', 'Save spike detection threshold', 'Callback', @(hO, hE) obj.saveThresh());

            % Edit menu
            editMenu = obj.hFigPreview.figApply(@uimenu, 'Label', 'Edit');
            uimenu(editMenu, 'Label', 'Bad site threshold', 'Callback', @(hO, hE) obj.setSiteCorrThresh());
            uimenu(editMenu, 'Label', 'Spike detection threshold', 'Callback', @(hO, hE) obj.setQQFactor());

            filterMenu = uimenu(editMenu, 'Label', 'Filter mode');
            obj.menuOptions(filterMenu, {'ndiff', 'bandpass', 'sgdiff', 'fir1', 'user'}, @obj.setFilter);
            obj.menuCheckbox(filterMenu, obj.hCfg.filterType);

            refMenu = uimenu(editMenu, 'Label', 'Reference mode');
            obj.menuOptions(refMenu, {'none', 'mean', 'median'}, @obj.setCARMode); % @TODO: local mean
            uimenu(editMenu, 'Label', 'Common reference threshold', 'Callback', @(hO, hE) obj.setBlankThresh());
            uimenu(editMenu, 'Label', 'FFT cleanup threshold', 'Callback', @(hO, hE) obj.setFFTThreshMAD());

            % View menu
            viewMenu = obj.hFigPreview.figApply(@uimenu, 'Label', 'View');

            tRangeMenu = uimenu(viewMenu, 'Label', 'Display time range (s)');
            obj.menuOptions(tRangeMenu, {'0.05', '0.1', '0.2', '0.5', '1', '2', '5', 'Custom'}, @obj.setTimeRange);
            obj.menuCheckbox(tRangeMenu, '0.1');

            uimenu(viewMenu, 'Label', 'Display site range', 'Callback', @(hO, hE) obj.setSiteRange());
            uimenu(viewMenu, 'Label', 'Show raw traces [F]', 'Callback', @(hO, hE) obj.keyPressFigPreview([], struct('Key', 'f')), 'Tag', 'menu_preview_view_filter');
            uimenu(viewMenu, 'Label', 'Show spike [T]hreshold', 'Callback', @(hO, hE) obj.keyPressFigPreview([], struct('Key', 't')), 'Tag', 'menu_preview_view_threshold');
            uimenu(viewMenu, 'Label', 'Toggle [S]pikes', 'Callback', @(hO, hE) obj.keyPressFigPreview([], struct('Key', 's')), 'Tag', 'menu_preview_view_spike');
            uimenu(viewMenu, 'Label', 'Toggle [G]rid', 'Callback', @(hO, hE) obj.keyPressFigPreview([], struct('Key', 'g')), 'Tag', 'menu_preview_view_grid');

            siteMenu = uimenu(viewMenu, 'Label', 'Site view');
            obj.menuOptions(siteMenu, {'Site correlation', 'Spike threshold', 'Event rate (Hz)', 'Event SNR (median)'}, @obj.setSiteView);
            obj.menuCheckbox(siteMenu, 'Site correlation');

            refMenu = uimenu(viewMenu, 'Label', 'Reference view');
            obj.menuOptions(refMenu, {'original', 'binned'}, @obj.setRefView);
            obj.menuCheckbox(refMenu, 'original');

            freqMenu = uimenu(viewMenu, 'Label', 'Frequency scale');
            obj.menuOptions(freqMenu, {'Linear', 'Log'}, @obj.setFreqScale);
            obj.menuCheckbox(freqMenu, 'Linear');

            % mh_view_psd = uimenu(viewMenu, 'Label', 'PSD view');
            % obj.menuOptions(mh_view_psd, {'original', 'detrended'}, @Fig_preview_view_psd_);
            % obj.menuCheckbox(mh_view_psd, 'original');
        end

        function exportTraces(obj)
            %EXPORTTRACES Export raw and filtered traces
            jrclust.utils.exportToWorkspace(struct('tracesFilt', obj.tracesFilt, ...
                                                   'tracesRaw', obj.tracesRaw), 1);
        end

        function loadPreview(obj)
            %LOADPREVIEW
            nLoadsMax = obj.hCfg.nLoadsMaxPreview;
            nSecsLoad = obj.hCfg.nSecsLoadPreview;

            % load files
            rawRecordings = jrclust.utils.subsample(obj.hCfg.rawRecordings, nLoadsMax);
            nFiles = numel(rawRecordings);
            obj.nLoadsPerFile = floor(nLoadsMax / nFiles);
            obj.nSamplesPerLoad = ceil(nSecsLoad*obj.hCfg.sampleRate);

            obj.hCfg.useGPU = 0;

            tracesRaw_ = cell(nFiles, 1);

            for iFile = 1:nFiles
                hRec = jrclust.detect.newRecording(rawRecordings{iFile}, obj.hCfg);
                if hRec.nSamples <= obj.nSamplesPerLoad
                    nLoadsFile = 1;
                    nSamplesLoad = hRec.nSamples;
                else
                    nLoadsFile = min(obj.nLoadsPerFile, floor(hRec.nSamples / obj.nSamplesPerLoad));
                    nSamplesLoad = obj.nSamplesPerLoad;
                end

                multiBounds = jrclust.views.sampleSkip([1, nSamplesLoad], hRec.nSamples, nLoadsFile);
                obj.fileBounds(hRec.rawPath) = multiBounds;

                fileTraces_ = cell(nLoadsFile, 1);

                hRec.openRaw();
                for iLoad = 1:nLoadsFile
                    iBounds = multiBounds{iLoad};
                    fileTraces_{iLoad} = hRec.readRawROI(obj.hCfg.siteMap, iBounds(1):iBounds(2))';
                end
                hRec.closeRaw();

                tracesRaw_{iFile} = cat(1, fileTraces_{:});
            end

            obj.tracesRaw = cat(1, tracesRaw_{:});
        end

        function saveConfig(obj)
            %SAVEPRM Update parameter file, show a preview window
            if isempty(obj.figData)
                return;
            end

            % Select variables to export
            newPrms = struct('fftThresh', obj.fftThresh, ...
                             'ignoreSites', find(obj.ignoreMe), ...
                             'qqFactor', obj.qqFactor, ...
                             'filterType', obj.filterType, ...
                             'CARMode', obj.CARMode, ...
                             'blankThresh', obj.blankThresh, ...
                             'blankPeriod', obj.blankPeriod);

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

            obj.hCfg.fftThresh = newPrms.fftThresh;
            obj.hCfg.ignoreSites = newPrms.ignoreSites;
            obj.hCfg.qqFactor = newPrms.qqFactor;
            obj.hCfg.filterType = newPrms.filterType;
            obj.hCfg.CARMode = newPrms.CARMode;
            obj.hCfg.blankThresh = newPrms.blankThresh;
            obj.hCfg.blankPeriod = newPrms.blankPeriod;

            obj.hCfg.save();
            obj.hCfg.edit();
        end

        function saveThresh(obj)
            %SAVETHRESH Export siteThresh to file, update hCfg.threshFile
            siteThresh = obj.siteThresh; %#ok<NASGU,PROP> (for saving)
            threshFile = jrclust.utils.subsExt(obj.hCfg.configFile, '_thresh.mat');
            save(threshFile, 'siteThresh'); % also need to store filter values?

            obj.hCfg.threshFile = threshFile;
            obj.hCfg.save();
            obj.hCfg.edit();

            jrclust.utils.qMsgBox(sprintf('Saved to %s and updated %s (threshFile)', ...
                threshFile, obj.hCfg.configFile));
        end

        function setBlankThresh(obj)
            %SETBLANKTHRESH
            blankThresh_ = num2str(obj.blankThresh);
            blank_period_ms_ = num2str(obj.blankPeriod);

            dlgAns = inputdlg({'blankThresh (MAD)', 'blankPeriod (millisecond)'}, 'Common reference threshold', 1, {blankThresh_, blank_period_ms_});
            if isempty(dlgAns)
                return;
            end

            blankThresh_ = str2double(dlgAns{1});
            blank_period_ms_ = str2double(dlgAns{2});

            if isempty(blankThresh_)
                blankThresh_ = 0;
            end
            if isnan(blankThresh_) || isnan(blank_period_ms_)
                return;
            end

            obj.blankThresh = blankThresh_;
            obj.blankPeriod = blank_period_ms_;

            obj.updateFigPreview(1);
        end

        function setCARMode(obj, label, hMenu)
            %SETCARMODE Update CAR mode
            obj.menuCheckbox(hMenu, label);
            obj.CARMode = label;

            obj.updateFigPreview(1);
        end

        function setFFTThreshMAD(obj)
            %SETFFTTHRESHMAD
            if isempty(obj.fftThresh)
                fftThresh_ = '0';
            else
                fftThresh_ = num2str(obj.fftThresh);
            end

            dlgAns = inputdlg('Set a threshold for FFT cleanup (0 to disable)', 'fftThresh (20 recommended)', 1, {fftThresh_});
            if isempty(dlgAns)
                return;
            end

            fftThresh_ = str2double(dlgAns{1});
            if isnan(fftThresh_)
                return;
            end

            obj.fftThresh = fftThresh_;
            obj.updateFigPreview(1);
        end

        function setFilter(obj, label, hMenu)
            %SETFILTER Update filter
            obj.menuCheckbox(hMenu, label);
            obj.filterType = label;
            obj.fFilter = 1;

            obj.updateFigPreview(1);
        end

        function setFreqScale(obj, label, hMenu)
            %SETFREQSCALE
            obj.menuCheckbox(hMenu, label);
            switch label
                % case 'Power'
                % case 'Detrended'
                case 'Linear'
                    obj.hFigPreview.axApply('hAxPSD', @set, 'XScale', 'linear');

                case 'Log'
                    obj.hFigPreview.axApply('hAxPSD', @set, 'XScale', 'log');
            end
        end

        function setQQFactor(obj)
            % update spike threshold qqFactor
            dlgAns = inputdlg('Set a spike detection threshold (qqFactor)', 'Spike detection threshold', 1, {num2str(obj.qqFactor)});
            if isempty(dlgAns)
                return;
            end

            qqFactor_ = str2double(dlgAns{1});
            if isnan(qqFactor_) || isempty(qqFactor_)
                return;
            end

            obj.qqFactor = qqFactor_;
            obj.updateFigPreview(1);
        end

        function setRefView(obj, label, hMenu)
            %SETREFVIEW
            obj.refView = label;
            obj.menuCheckbox(hMenu, label);
            doPlotFigPreview(obj.hFigPreview, obj.figData, 1, obj.hCfg);
        end

        function setSiteCorrThresh(obj)
            %SETSITECORRTHRESH Update the siteCorr threshold to determine a bad site
            if isempty(obj.siteCorrThresh)
                siteCorrThreshStr = '0';
            else
                siteCorrThreshStr = num2str(obj.siteCorrThresh);
            end

            dlgAns = inputdlg('Set a correlation threshold [0-1) for detecting bad sites (0 to disable)', 'siteCorrThresh (set to 0 to disable)', 1, {siteCorrThreshStr});
            if isempty(dlgAns)
                return;
            end

            siteCorrThresh_ = str2double(dlgAns{1});
            if siteCorrThresh_ >= 1 || siteCorrThresh_ < 0 || isnan(siteCorrThresh_)
                jrclust.utils.qMsgBox(sprintf('Invalid range: %s', dlgAns{1}));
                return;
            end

            obj.siteCorrThresh = siteCorrThresh_;
            obj.updateFigPreview(1);
        end

        function setSiteRange(obj)
            %SETSITERANGE
            siteLo = num2str(obj.siteLim(1));
            siteHi = num2str(obj.siteLim(end));

            dlgAns = inputdlg({'Show site from (>=1)', sprintf('Show site to (<=%d)', obj.hCfg.nSites)}, 'Display site range', 1, {siteLo, siteHi});
            if isempty(dlgAns)
                return;
            end

            siteStart = max(str2double(dlgAns{1}), 1);
            siteEnd = min(str2double(dlgAns{2}), obj.hCfg.nSites);
            obj.siteLim = [siteStart, siteEnd];

            obj.hFigPreview.axApply('hAxTraces', @set, 'YLim', obj.siteLim + [-1, 1]);
            obj.hFigPreview.axApply('hAxSites', @set, 'YLim', obj.siteLim + [-1, 1]);
        end

        function setSiteView(obj, label, hMenu)
            %SETSITEVIEW
            obj.siteView = label;
            obj.menuCheckbox(hMenu, label);
            doPlotFigPreview(obj.hFigPreview, obj.figData, 1, obj.hCfg);
        end

        function setTimeRange(obj, label, hMenu)
            %SETTIMERANGE Sets a display time range
            if nargin < 2
                label = 'Custom';
            end

            if strcmp(label, 'Custom') % ask user input box
                dlgAns = inputdlg('Display time range (s)', 'Time range in seconds', 1, {'.2'});
                if isempty(dlgAns)
                    return;
                end

                tRange = str2double(dlgAns{1});
            else
                tRange = str2double(label);
            end

            if isnan(tRange)
                return;
            end
            obj.menuCheckbox(hMenu, label);

            obj.windowWidth = round(tRange * obj.hCfg.sampleRate);
            obj.windowBounds = obj.windowBounds(1) + [0, obj.windowWidth - 1];
            doPlotFigPreview(obj.hFigPreview, obj.figData, 0, obj.hCfg);
        end

        function updateFigPreview(obj, fKeepView)
            %UPDATEFIGPREVIEW Update parameters and replot
            obj.hFigPreview.wait(1);

            if obj.fftThresh > 0
                obj.tracesClean = jrclust.filters.fftClean(obj.tracesRaw, obj.fftThresh, obj.hCfg); % fft filter
            else
                obj.tracesClean = obj.tracesRaw;
            end

            % Find bad sites
            if obj.siteCorrThresh > 0
                obj.ignoreMe = (obj.maxCorrSite < obj.siteCorrThresh);
            elseif ~isempty(obj.hCfg.ignoreSites)
                obj.ignoreMe = false(size(obj.maxCorrSite));
                obj.ignoreMe(obj.hCfg.ignoreSites) = 1;
            else
                obj.ignoreMe = false(size(obj.maxCorrSite));
            end

            % filter and CAR
            obj.hCfg.setTemporaryParams('CARMode', 'none', 'useGPU', 0, 'filterType', obj.filterType, ...
                'blankPeriod', obj.blankPeriod, 'blankThresh', obj.blankThresh, 'useParfor', 0);
            obj.tracesFilt = jrclust.filters.filtCAR(obj.tracesClean, [], [], 0, obj.hCfg);
            obj.tracesCAR = jrclust.utils.getCAR(obj.tracesFilt, obj.CARMode, obj.ignoreSites);

            if ~strcmpi(obj.CARMode, 'none')
                obj.tracesFilt = bsxfun(@minus, obj.tracesFilt, cast(obj.tracesCAR, 'like', obj.tracesFilt));
            end

            obj.tracesCAR = jrclust.utils.madScore(mean(obj.tracesCAR, 2)); % Save in MAD unit

            % reference threshold
            [rawPSD, obj.psdFreq] = getPSD(obj.tracesRaw(:, ~obj.ignoreMe), obj.hCfg.sampleRate, 4);
            cleanPSD = getPSD(obj.tracesClean(:, ~obj.ignoreMe), obj.hCfg.sampleRate, 4);

            obj.psdPower = mean(rawPSD, 2);
            obj.psdPowerClean = mean(cleanPSD, 2);

            % find threshold and detect spikes
            siteRMS = jrclust.utils.estimateRMS(obj.tracesFilt, 1e5);
            obj.siteThresh = int16(siteRMS * obj.qqFactor);
            obj.siteThresh(obj.ignoreMe) = 0;

            obj.isThreshCrossing = bsxfun(@lt, obj.tracesFilt, -abs(obj.siteThresh));
            obj.isThreshCrossing(:, obj.ignoreMe) = 0; % ignore threshold crossings on bad sites

            % Spike detection
            [obj.keepMe, obj.channelMeansMAD] = jrclust.utils.carReject(obj.tracesCAR, obj.hCfg.blankPeriod, obj.hCfg.blankThresh, obj.hCfg.sampleRate);
            [obj.spikeTimes, obj.spikeAmps, obj.spikeSites] = jrclust.utils.detectPeaks(obj.tracesFilt, obj.siteThresh, obj.keepMe, obj.hCfg);

            durationSecs = size(obj.tracesFilt, 1) / obj.hCfg.sampleRate;

            obj.siteEventRate = hist(obj.spikeSites, 1:obj.hCfg.nSites)/durationSecs; % event count
            obj.siteEventSNR = abs(single(arrayfun(@(i) median(obj.spikeAmps(obj.spikeSites == i)), 1:obj.hCfg.nSites)))./siteRMS;

            obj.hCfg.resetTemporaryParams();

            doPlotFigPreview(obj.hFigPreview, obj.figData, fKeepView, obj.hCfg);
            obj.hFigPreview.wait(0);
        end
    end

    %% USER METHODS
    methods
        function preview(obj)
            %PREVIEW show preview plot
            obj.loadPreview();

            nSites = size(obj.tracesRaw, 2);

            siteCorr = corr(single(obj.tracesRaw));
            siteCorr(logical(eye(size(siteCorr)))) = 0;

            % Create a Figure
            gap = .05;
            obj.hFigPreview = jrclust.views.Figure('FigPreview', [0 0 .5 1], obj.hCfg.configFile, 1, 1); %plot a summary pannel
            obj.hFigPreview.addAxes('hAxMean', 'Position', [gap, gap, 3/4-gap, 1/4-gap], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxTraces', 'Position', [gap, 1/4+gap, 3/4-gap, 3/4-gap*2], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxSites', 'Position', [3/4+gap, gap, 1/4-gap*1.5, 2/3-gap*2], 'NextPlot', 'add');
            obj.hFigPreview.addAxes('hAxPSD', 'Position', [3/4+gap, 2/3+gap, 1/4-gap*1.5, 1/3-gap*2], 'NextPlot', 'add');
            linkaxes([obj.hFigPreview.hAxes('hAxMean'), obj.hFigPreview.hAxes('hAxTraces')], 'x');

            % Callback functions
            obj.addMenu();
            obj.hFigPreview.figApply(@set, 'KeyPressFcn', @obj.keyPressFigPreview, 'BusyAction', 'cancel');
            %mouse_figure(obj.hFigPreview, hAxTraces);
            obj.hFigPreview.setMouseable([], 'hAxTraces');

            % Build figData
            obj.filterType = obj.hCfg.filterType;
            obj.CARMode = obj.hCfg.CARMode;
            obj.siteCorrThresh = obj.hCfg.siteCorrThresh;
            obj.fftThresh = obj.hCfg.fftThresh;
            obj.qqFactor = obj.hCfg.qqFactor;
            obj.blankThresh = obj.hCfg.blankThresh;
            obj.blankPeriod = obj.hCfg.blankPeriod;

            ignoreMe_ = false(size(obj.hCfg.siteMap));
            ignoreMe_(obj.hCfg.ignoreSites) = 1;
            obj.ignoreMe = ignoreMe_;

            obj.maxAmp = obj.hCfg.maxAmp;
            obj.maxCorrSite = max(siteCorr);
            obj.siteLim = [1 nSites];
            obj.fFilter = 1;
            obj.fGrid = 1;
            obj.fShowThresh = 0;
            obj.fShowSpikes = 1;
            obj.siteView = 'Site correlation';
            obj.refView = 'binned';
            obj.psdView = 'original';

            if isprop(obj.hCfg, 'preview_window')
                obj.windowWidth = round(obj.hCfg.preview_window * obj.hCfg.sampleRate); % TW
            else
                obj.windowWidth = obj.nSamplesPerLoad;
            end

            obj.windowBounds = [1, obj.windowWidth];
            obj.helpText = {'Left/Right: change time (Shift: x4)', ...
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

            drawnow;
            obj.updateFigPreview(0);
        end
    end

    %% GETTERS/SETTERS
    methods
        % blankThresh
        function bt = get.blankThresh(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'blankThresh')
                bt = [];
            else
                bt = obj.figData.blankThresh;
            end
        end
        function set.blankThresh(obj, bt)
            obj.figData.blankThresh = bt;
        end

        % blankPeriod
        function bp = get.blankPeriod(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'blankPeriod')
                bp = [];
            else
                bp = obj.figData.blankPeriod;
            end
        end
        function set.blankPeriod(obj, bp)
            obj.figData.blankPeriod = bp;
        end

        % CARMode
        function cm = get.CARMode(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'CARMode')
                cm = [];
            else
                cm = obj.figData.CARMode;
            end
        end
        function set.CARMode(obj, cm)
            obj.figData.CARMode = cm;
        end

        % channelMeansMAD
        function cm = get.channelMeansMAD(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'channelMeansMAD')
                cm = [];
            else
                cm = obj.figData.channelMeansMAD;
            end
        end
        function set.channelMeansMAD(obj, cm)
            obj.figData.channelMeansMAD = cm;
        end

        % fFilter
        function ff = get.fFilter(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'fFilter')
                ff = [];
            else
                ff = obj.figData.fFilter;
            end
        end
        function set.fFilter(obj, ff)
            obj.figData.fFilter = ff;
        end

        % fftThresh
        function ft = get.fftThresh(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'fftThresh')
                ft = [];
            else
                ft = obj.figData.fftThresh;
            end
        end
        function set.fftThresh(obj, ft)
            obj.figData.fftThresh = ft;
        end

        % fGrid
        function fg = get.fGrid(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'fGrid')
                fg = [];
            else
                fg = obj.figData.fGrid;
            end
        end
        function set.fGrid(obj, fg)
            obj.figData.fGrid = fg;
        end

        % filterType
        function ft = get.filterType(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'filterType')
                ft = [];
            else
                ft = obj.figData.filterType;
            end
        end
        function set.filterType(obj, ft)
            obj.figData.filterType = ft;
        end

        % fShowSpikes
        function fs = get.fShowSpikes(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'fShowSpikes')
                fs = [];
            else
                fs = obj.figData.fShowSpikes;
            end
        end
        function set.fShowSpikes(obj, fs)
            obj.figData.fShowSpikes = fs;
        end

        % fShowThresh
        function fs = get.fShowThresh(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'fShowThresh')
                fs = [];
            else
                fs = obj.figData.fShowThresh;
            end
        end
        function set.fShowThresh(obj, fs)
            obj.figData.fShowThresh = fs;
        end

        % helpText
        function ht = get.helpText(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'helpText')
                ht = [];
            else
                ht = obj.figData.helpText;
            end
        end
        function set.helpText(obj, ht)
            obj.figData.helpText = ht;
        end

        % ignoreMe
        function ig = get.ignoreMe(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'ignoreMe')
                ig = [];
            else
                ig = obj.figData.ignoreMe;
            end
        end
        function set.ignoreMe(obj, ig)
            obj.figData.ignoreMe = ig;
            obj.figData.ignoreSites = find(ig);
        end

        % ignoreSites
        function is = get.ignoreSites(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'ignoreSites')
                is = [];
            else
                is = obj.figData.ignoreMe;
            end
        end

        % isThreshCrossing
        function it = get.isThreshCrossing(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'isThreshCrossing')
                it = [];
            else
                it = obj.figData.isThreshCrossing;
            end
        end
        function set.isThreshCrossing(obj, it)
            obj.figData.isThreshCrossing = it;
        end

        % keepMe
        function km = get.keepMe(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'keepMe')
                km = [];
            else
                km = obj.figData.keepMe;
            end
        end
        function set.keepMe(obj, km)
            obj.figData.keepMe = km;
        end

        % maxAmp
        function ma = get.maxAmp(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'maxAmp')
                ma = [];
            else
                ma = obj.figData.maxAmp;
            end
        end
        function set.maxAmp(obj, ma)
            obj.figData.maxAmp = ma;
        end

        % maxCorrSite
        function mc = get.maxCorrSite(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'maxCorrSite')
                mc = [];
            else
                mc = obj.figData.maxCorrSite;
            end
        end
        function set.maxCorrSite(obj, mc)
            obj.figData.maxCorrSite = mc;
        end

        % nLoads
        function nl = get.nLoads(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'nLoads')
                nl = [];
            else
                nl = obj.figData.nLoads;
            end
        end
        function set.nLoads(obj, nl)
            obj.figData.nLoads = nl;
        end

        % nSamplesTotal
        function ns = get.nSamplesTotal(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'tracesRaw')
                ns = [];
            else
                ns = size(obj.tracesRaw, 1);
            end
        end

        % psdFreq
        function pf = get.psdFreq(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'psdFreq')
                pf = [];
            else
                pf = obj.figData.psdFreq;
            end
        end
        function set.psdFreq(obj, pf)
            obj.figData.psdFreq = pf;
        end

        % psdPower
        function pp = get.psdPower(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'psdPower')
                pp = [];
            else
                pp = obj.figData.psdPower;
            end
        end
        function set.psdPower(obj, pp)
            obj.figData.psdPower = pp;
        end

        % psdPowerClean
        function pp = get.psdPowerClean(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'psdPowerClean')
                pp = [];
            else
                pp = obj.figData.psdPowerClean;
            end
        end
        function set.psdPowerClean(obj, pp)
            obj.figData.psdPowerClean = pp;
        end

        % psdView
        function pv = get.psdView(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'psdView')
                pv = [];
            else
                pv = obj.figData.psdView;
            end
        end
        function set.psdView(obj, pv)
            obj.figData.psdView = pv;
        end

        % qqFactor
        function qf = get.qqFactor(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'qqFactor')
                qf = [];
            else
                qf = obj.figData.qqFactor;
            end
        end
        function set.qqFactor(obj, pp)
            obj.figData.qqFactor = pp;
        end

        % refView
        function rv = get.refView(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'refView')
                rv = [];
            else
                rv = obj.figData.refView;
            end
        end
        function set.refView(obj, rv)
            obj.figData.refView = rv;
        end

        % siteCorrThresh
        function sc = get.siteCorrThresh(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteCorrThresh')
                sc = [];
            else
                sc = obj.figData.siteCorrThresh;
            end
        end

        function set.siteCorrThresh(obj, sc)
            obj.figData.siteCorrThresh = sc;
        end

        % siteEventRate
        function se = get.siteEventRate(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteEventRate')
                se = [];
            else
                se = obj.figData.siteEventRate;
            end
        end
        function set.siteEventRate(obj, se)
            obj.figData.siteEventRate = se;
        end

        % siteEventSNR
        function se = get.siteEventSNR(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteEventSNR')
                se = [];
            else
                se = obj.figData.siteEventSNR;
            end
        end
        function set.siteEventSNR(obj, se)
            obj.figData.siteEventSNR = se;
        end

        % siteSNR
        function ss = get.siteSNR(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteSNR')
                ss = [];
            else
                ss = obj.figData.siteSNR;
            end
        end
        function set.siteSNR(obj, ss)
            obj.figData.siteSNR = ss;
        end

        % siteThresh
        function st = get.siteThresh(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteThresh')
                st = [];
            else
                st = obj.figData.siteThresh;
            end
        end
        function set.siteThresh(obj, st)
            obj.figData.siteThresh = st;
        end

        % siteLim
        function sl = get.siteLim(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteLim')
                sl = [];
            else
                sl = obj.figData.siteLim;
            end
        end
        function set.siteLim(obj, sl)
            obj.figData.siteLim = max(1, min(obj.hCfg.nSites, sl));
        end

        % siteView
        function sv = get.siteView(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'siteView')
                sv = [];
            else
                sv = obj.figData.siteView;
            end
        end
        function set.siteView(obj, sv)
            obj.figData.siteView = sv;
        end

        % spikeAmps
        function sa = get.spikeAmps(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'spikeAmps')
                sa = [];
            else
                sa = obj.figData.spikeAmps;
            end
        end
        function set.spikeAmps(obj, sa)
            obj.figData.spikeAmps = sa;
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'spikeSites')
                ss = [];
            else
                ss = obj.figData.spikeSites;
            end
        end
        function set.spikeSites(obj, ss)
            obj.figData.spikeSites = ss;
        end

        % spikeTimes
        function st = get.spikeTimes(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'spikeTimes')
                st = [];
            else
                st = obj.figData.spikeTimes;
            end
        end
        function set.spikeTimes(obj, st)
            obj.figData.spikeTimes = st;
        end

        % tracesCAR
        function tc = get.tracesCAR(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'tracesCAR')
                tc = [];
            else
                tc = obj.figData.tracesCAR;
            end
        end
        function set.tracesCAR(obj, tc)
            obj.figData.tracesCAR = tc;
        end

        % tracesClean
        function tc = get.tracesClean(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'tracesClean')
                tc = [];
            else
                tc = obj.figData.tracesClean;
            end
        end
        function set.tracesClean(obj, tc)
            obj.figData.tracesClean = tc;
        end

        % tracesFilt
        function tc = get.tracesFilt(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'tracesFilt')
                tc = [];
            else
                tc = obj.figData.tracesFilt;
            end
        end
        function set.tracesFilt(obj, tc)
            obj.figData.tracesFilt = tc;
        end

        % tracesRaw
        function tc = get.tracesRaw(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'tracesRaw')
                tc = [];
            else
                tc = obj.figData.tracesRaw;
            end
        end
        function set.tracesRaw(obj, tc)
            obj.figData.tracesRaw = tc;
        end

        % windowBounds
        function wb = get.windowBounds(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'windowBounds')
                wb = [];
            else
                wb = obj.figData.windowBounds;
            end
        end
        function set.windowBounds(obj, wb)
            if wb(1) < 1
                wb = [1, obj.windowWidth];
            elseif wb(2) > obj.nSamplesTotal
                wb = [1 - obj.windowWidth, 0] + obj.nSamplesTotal;
            end

            obj.figData.windowBounds = wb;
        end

        % windowWidth
        function ww = get.windowWidth(obj)
            if isempty(obj.figData) || ~isfield(obj.figData, 'windowWidth')
                ww = [];
            else
                ww = obj.figData.windowWidth;
            end
        end
        function set.windowWidth(obj, ww)
            ww = min(ww, obj.nSamplesTotal);
            obj.figData.windowWidth = ww;
        end
    end
end
