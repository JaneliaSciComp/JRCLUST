classdef TracesController < handle
    %TRACESCONTROLLER Summary of this class goes here
    
    properties (SetAccess=private)
        hCfg;
        hClust;
        hRec;
    end

    properties (SetAccess=private, Transient)
        hFigPSD;
        hFigTraces;

        showFilter;
        showGrid;
        showSpikes;
        showTraces;

        tracesFull;
        tracesRaw;
        tracesFilt;
    end

    %% LIFECYCLE
    methods
        function obj = TracesController(hCfg)
            %TRACESCONTROLLER Construct an instance of this class
            obj.hCfg = hCfg;
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods
        function keyPressFigTraces(obj, hObject, hEvent)
            % 2017/6/22 James Jun: Added nTime_traces multiview
            factor = 4^double(keyMod(hEvent, 'shift')); % 1 or 4
            nSites = numel(obj.hCfg.siteMap);

            switch hEvent.Key
                case 'uparrow'
                    obj.hFigTraces.figData.maxAmp = obj.hFigTraces.figData.maxAmp*sqrt(2)^-factor;
                    obj.updateFigTraces(false);

                case 'downarrow'
                    obj.hFigTraces.figData.maxAmp = obj.hFigTraces.figData.maxAmp*sqrt(2)^factor;
                    obj.updateFigTraces(false);

                case 'e' %export current view
                    jrclust.utils.exportToWorkspace(struct('tracesRaw', obj.tracesRaw, 'tracesFilt', obj.tracesFilt), true);

                case 'f' % toggle filter
                    obj.toggleFilter();

                case 'g' % toggle grid
                    obj.toggleGrid();

                case 'h'
                    msgbox_(obj.hFigTraces.figData.helpText, 1);

                case 'p' % power spectrum
                    iSite = inputdlgNum(sprintf('Site# to show (1-%d, 0 for all)', nSites), 'Site#', 0);

                    if isnan(iSite)
                        return;
                    end

                    obj.hFigPSD = jrclust.views.Figure('FigPsd', [.5 0 .5 1], obj.hCfg.configFile, true, true); % show to the right

                    % ask user which channels to plot
                    if iSite > 0
                        tracesFilt_ = obj.tracesFilt(:, iSite)';
                    else
                        tracesFilt_ = obj.tracesFilt;
                    end

                    doPlotFigPSD(obj.hFigPSD, tracesFilt_, obj.hCfg); %obj.hCfg.sampleRate/obj.hCfg.nSkip_show, 'ignoreSites', obj.hCfg.ignoreSites));

                case 'r' % reset view
                    if obj.hFigTraces.figData.maxAmp ~= obj.hCfg.maxAmp
                        obj.hFigTraces.figData.maxAmp = obj.hCfg.maxAmp;
                        obj.updateFigTraces(true);
                    else
                        resetFigTraces(obj.hFigTraces, obj.tracesRaw, obj.hCfg);
                    end

                case 's' %show/hide spikes
                    obj.toggleSpikes();

                case 't' %show/hide traces
                    obj.toggleTraces();
            end % switch

%             switch lower(hEvent.Key)
%                 
% 
% 
%                 case {'leftarrow', 'rightarrow', 'j', 'home', 'end'}
%                 switch lower(hEvent.Key)
%                     case 'leftarrow'
%                     nlim_bin = obj.hFigTraces.figData.nlim_bin - (obj.hFigTraces.figData.nLoad_bin) * factor; %no overlap
%                     if nlim_bin(1)<1
%                         msgbox_('Beginning of file', 1);
%                         nlim_bin = [1, obj.hFigTraces.figData.nLoad_bin];
%                     end
%                     case 'rightarrow'
%                     nlim_bin = obj.hFigTraces.figData.nlim_bin + (obj.hFigTraces.figData.nLoad_bin + 1) * factor; %no overlap
%                     if nlim_bin(2) > obj.hFigTraces.figData.nSamples_bin
%                         msgbox_('End of file', 1);
%                         nlim_bin = [-obj.hFigTraces.figData.nLoad_bin+1, 0] + obj.hFigTraces.figData.nSamples_bin;
%                     end
%                     case 'home' %beginning of file
%                     nlim_bin = [1, obj.hFigTraces.figData.nLoad_bin];
%                     case 'end' %end of file
%                     nlim_bin = [-obj.hFigTraces.figData.nLoad_bin+1, 0] + obj.hFigTraces.figData.nSamples_bin;
%                     case 'j'
%                     vcAns = inputdlg('Go to time (s)', 'Jump to time', 1, {'0'});
%                     if isempty(vcAns), return; end
%                     try
%                         nlim_bin = round(str2double(vcAns)*obj.hCfg.sRateHz) + [1, obj.hFigTraces.figData.nLoad_bin];
%                     catch
%                         return;
%                     end
%                 end %switch
%                 nTime_traces = get_(obj.hCfg, 'nTime_traces');
%                 [cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, obj.hFigTraces.figData.nSamples_bin, nTime_traces);
%                 if obj.hCfg.fTranspose_bin
%                     fseek_(obj.hFigTraces.figData.fid_bin, nlim_bin(1), obj.hCfg);
%                     if nTime_traces > 1
%                         obj.tracesRaw = load_bin_multi_(obj.hFigTraces.figData.fid_bin, cvn_lim_bin, obj.hCfg)';
%                     else
%                         obj.tracesRaw = jrclust.utils.readBin(obj.hFigTraces.figData.fid_bin, obj.hCfg.vcDataType, [obj.hCfg.nChans, obj.hFigTraces.figData.nLoad_bin])';
%                     end
%                 else
%                     obj.tracesRaw = obj.tracesFull(viRange_bin, :);
%                 end
%                 obj.tracesRaw = u2i(obj.tracesRaw);
%                 obj.hFigTraces.figData.nlim_bin = nlim_bin;
%                 doPlotFigTraces(1); %redraw
% 
% 
%                 case 'c' %channel query
%                 msgbox_('Draw a rectangle', 1);
%                 obj.hRect = imrect_();
%                 if isempty(obj.hRect), return ;end
%                 vrPos_rect = getPosition(obj.hRect);
%                 S_plot = get(obj.hFigTraces.figData.hPlot, 'UserData');
%                 vrX = get(obj.hFigTraces.figData.hPlot, 'XData');
%                 vrY = get(obj.hFigTraces.figData.hPlot, 'YData');
%                 viIndex = find(vrX >= vrPos_rect(1) & vrX <= sum(vrPos_rect([1,3])) & vrY >= vrPos_rect(2) & vrY <= sum(vrPos_rect([2,4])));
%                 if isempty(viIndex), delete_multi_(obj.hRect); return; end
%                 index_plot = round(median(viIndex));
%                 [~, iSite] = ind2sub(size(obj.tracesFilt), index_plot);
%                 mrX = reshape(vrX, S_plot.dimm);
%                 mrY = reshape(vrY, S_plot.dimm);
%                 hold(obj.hFigTraces.figData.hAx, 'on');
%                 hPoint = plot(vrX(index_plot), vrY(index_plot), 'r*');
%                 hLine = plot(obj.hFigTraces.figData.hAx, mrX(:,iSite), mrY(:,iSite), 'r-');
%                 hold(obj.hFigTraces.figData.hAx, 'off');
%                 iChan = obj.hCfg.viSite2Chan(iSite);
%                 msgbox_(sprintf('Site: %d/ Chan: %d', iSite, iChan), 1);
%                 delete_multi_(obj.hRect, hLine, hPoint);
%             end
        end
    end

    %% UTILITY METHODS
    methods (Hidden)
        function closeFigTraces(obj, hFig, hEvent)
            obj.hRec.close();
            delete(hFig);
        end

        function toggleFilter(obj)
            obj.showFilter = ~obj.showFilter;
            obj.hFigTraces.figData.filter = jrclust.utils.ifEq(obj.showFilter, 'on', 'off');
            obj.updateFigTraces(false);
        end

        function toggleGrid(obj)
            obj.showGrid = ~obj.showGrid;
            obj.hFigTraces.figData.grid = jrclust.utils.ifEq(obj.showGrid, 'on', 'off');
            obj.hFigTraces.axApply(@grid, obj.hFigTraces.figData.grid);
        end

        function toggleSpikes(obj)
            obj.showSpikes = ~obj.showSpikes;
            obj.hFigTraces.figData.spikes = jrclust.utils.ifEq(obj.showSpikes, 'on', 'off');
            obj.updateFigTraces(false);
        end

        function toggleTraces(obj)
            obj.showTraces = ~obj.showTraces;
            obj.hFigTraces.figData.traces = jrclust.utils.ifEq(obj.showTraces, 'on', 'off');
            obj.updateFigTraces(false);
        end

        function updateFigTraces(obj, resetAxes)
            doPlotFigTraces(obj.hFigTraces, obj.hCfg, obj.tracesRaw, resetAxes, obj.hClust);
        end
    end

    %% USER METHODS
    methods
        function show(obj, recID, showLFP, hClust)
            if nargin < 4
                hClust = [];
            end
            if nargin < 3
                showLFP = false;
            end
            if nargin < 2
                recID = '';
            end

            obj.hClust = hClust;

            % get file to show
            if numel(obj.hCfg.rawRecordings) > 1
                if isempty(recID)
                    arrayfun(@(i) fprintf('%d: %s\n', i, obj.hCfg.rawRecordings{i}), 1:numel(obj.hCfg.rawRecordings), 'UniformOutput', 0);
                    fprintf('---------------------------------------------\n');
                    recID = input('Please specify File ID from the list above:', 's');
                end

                if isempty(recID)
                    return;
                end

                iFile = str2double(recID);

                try
                    recFilename = obj.hCfg.rawRecordings{iFile};
                catch
                    return;
                end
            else
                recFilename = obj.hCfg.rawRecordings{1}; % if multiple files exist, load first
            end

            fprintf('Opening %s\n', recFilename);
            obj.hRec = jrclust.models.recording.Recording(recFilename, obj.hCfg);

        %     [fid_bin, nBytes_bin] = fopen_(vcFile_bin, 'r');
        %     if isempty(fid_bin)
        %         fprintf(2, '.bin file does not exist: %s\n', vcFile_bin);
        %         return;
        %     end
        %     nSamples_bin = floor(nBytes_bin / bytesPerSample_(obj.hCfg.vcDataType) / obj.hCfg.nChans);
            nSamplesTotal = obj.hRec.nSamples;
            windowWidth = min(floor(diff(obj.hCfg.dispTimeLimits) * obj.hCfg.sampleRate), nSamplesTotal);

            if obj.hCfg.dispTimeLimits(1) > 0
                iSampleOffset = ceil(obj.hCfg.dispTimeLimits(1) * obj.hCfg.sampleRate) + 1; % offset sample number
            else
                iSampleOffset = 1; %sample start location
            end

            windowBounds = [0, windowWidth - 1] + iSampleOffset;
            if windowBounds(1) < 1
                windowBounds = [1, windowWidth];
            end
            if windowBounds(2) > nSamplesTotal
                windowBounds = [-windowWidth + 1, 0] + nSamplesTotal;
            end

            multiBounds = sample_skip_(windowBounds, nSamplesTotal, obj.hCfg.nTime_traces);

            obj.tracesFull = [];
            tracesRaw_ = cellfun(@(lims) obj.hRec.readROI(obj.hCfg.siteMap, lims(1):lims(2)), multiBounds, 'UniformOutput', false);
            obj.tracesRaw = jrclust.utils.neCell2mat(tracesRaw_);

        %     if obj.hCfg.fTranspose_bin
        %         obj.tracesFull = [];
        %         fseek_(fid_bin, iSampleOffset, obj.hCfg);
        %         if obj.hCfg.nTime_traces > 1
        %             obj.tracesRaw = load_bin_multi_(fid_bin, cvn_lim_bin, obj.hCfg)';
        %         else
        %             obj.tracesRaw = jrclust.utils.readBin(fid_bin, obj.hCfg.vcDataType, [obj.hCfg.nChans, windowWidth])'; %next keypress: update tlim_show
        %         end
        %         %     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
        %     else %load whole thing
        %         obj.tracesFull = jrclust.utils.readBin(fid_bin, obj.hCfg.vcDataType, [nSamplesTotal, obj.hCfg.nChans]); %next keypress: update tlim_show
        %         fclose(fid_bin);
        %         fid_bin = [];
        %         %obj.tracesRaw = obj.tracesFull((nlim_bin(1):nlim_bin(2)), :);
        %         obj.tracesRaw = obj.tracesFull(viRange_bin, :);
        %         disp('Entire raw traces are cached to RAM since fTranspose=0.');
        %     end %if
            obj.tracesRaw = u2i(obj.tracesRaw);

            obj.hFigTraces = jrclust.views.Figure('FigTraces', [0 0 .5 1], recFilename, false, true);
            obj.hFigTraces.axes();
            obj.hFigTraces.addPlot('hLine', @line, nan, nan, 'Color', [1 1 1]*.5, 'LineWidth', .5);
            obj.hFigTraces.addPlot('hEdges', nan, nan, 'Color', [1 0 0]*.5, 'LineWidth', 1);
            obj.hFigTraces.axApply(@set, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');

            figData = struct('windowBounds', windowBounds, ...
                             'hRec', obj.hRec, ...
                             'nSamplesTotal', nSamplesTotal, ...
                             'windowWidth', windowWidth);

            figData.maxAmp = obj.hCfg.maxAmp;
            figData.title = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [C]han. query; [R]eset view; [P]SD; [S]pike; [A]ux chan; [E]xport; [T]race; [G]rid';
            figData.helpText = { ...
                'Left/Right: change time (Shift: x4)', ...
                '[J]ump T', ...
                '[Home/End]: go to beginning/end of file', ...
                '---------', ...
                'Up/Down: change scale (Shift: x4)', ...
                'Zoom: Mouse wheel', ...
                '[x/y/ESC]: zoom direction', ...
                'Pan: hold down the wheel and drag', ...
                '[R]eset view', ...
                '---------', ...
                '[F]ilter toggle', ...
                '[S]pike toggle', ...
                'Gri[D] toggle', ...
                '[T]races toggle', ...
                '---------', ...
                '[C]hannel query', ...
                '[A]ux channel display', ...
                '[P]ower spectrum', ...
                '[E]xport to workspace', ...
            };

            obj.showFilter = false;
            figData.filter = 'off';

            obj.showGrid = true;
            figData.grid = 'on';

            obj.showSpikes = false;
            figData.spikes = 'off';

            obj.showTraces = true;
            figData.traces = 'on';

            obj.hFigTraces.figData = figData;

            obj.hFigTraces.figApply(@set, 'Color', 'w', 'BusyAction', 'cancel', 'CloseRequestFcn', @obj.closeFigTraces);
            obj.hFigTraces.hFunKey = @obj.keyPressFigTraces;
            obj.hFigTraces.setMouseable();

            obj.tracesFilt = doPlotFigTraces(obj.hFigTraces, obj.hCfg, obj.tracesRaw, true, obj.hClust);
        end
    end
end

