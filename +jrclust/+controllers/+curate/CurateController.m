classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters

    properties (SetAccess=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hClust;         % Clustering object
    end

    properties (Dependent)
        hCfg;
    end

    properties (AbortSet, SetAccess=private, Hidden, Transient, SetObservable)
        currentSite;    % current site in FigTime
        hFigs;          % containers.Map of Figure objects
        hMenus;         % containers.Map of Menu handles
        isEnding;       % to prevent double prompting when endSession is called
        maxAmp;         % 
        projSites;      % current sites in FigProj
        selected;       % selected clusters, in order of selection
    end

    %% LIFECYCLE
    methods
        function obj = CurateController(hClust)
            obj.hClust = hClust;
            obj.hFigs = containers.Map();
            obj.hMenus = containers.Map();
            obj.currentSite = [];
            obj.maxAmp = [];
            obj.selected = [];
        end

        function delete(obj)
            if ~isempty(obj.hFigs)
                obj.closeFigures();
            end
        end
    end

    %% KEYPRESS/MOUSECLICK METHODS
    methods (Hidden)
        function keyPressFigSim(obj, hObject, hEvent) %#ok<*INUSL>
            %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
            hFigSim = obj.hFigs('FigSim');

            switch hEvent.Key
                case {'d', 'backspace', 'delete'} % delete
                    hFigSim.wait(true);
                    obj.deleteClusters();
                    hFigSim.wait(false);

                case 'k' % kilosort template similarity view
                    if obj.hCfg.fImportKsort && strcmp(hFigSim.figData.figView, 'waveform')
                        hFigSim.figData.figView = 'kilosort';
                        obj.updateFigSim();
                    end

                case 'm' % merge
                    hFigSim.wait(true);
                    obj.mergeSelected();
                    hFigSim.wait(false);

                case 's' % split
                    hFigSim.wait(true);
                    obj.autoSplit(true);
                    hFigSim.wait(false);

                case 'w' % waveform-based sim view
                    if obj.hCfg.fImportKsort && strcmp(hFigSim.figData.figView, 'kilosort')
                        hFigSim.figData.figView = 'waveform';
                        obj.updateFigSim();
                    end
            end % switch
        end

        function keyPressFigProj(obj, hObject, hEvent)
            %KEYPRESSFIGPROJ Handle callbacks for keys pressed in feature view
            hFigProj = obj.hFigs('FigProj');

            switch lower(hEvent.Key)
                case 'uparrow'
                    pow = -(4^double(keyMod(hEvent, 'shift'))); % 1 or 4
                    projScale = hFigProj.figData.boundScale*sqrt(2)^pow;
                    rescaleFigProj(hFigProj, projScale, obj.hCfg);

                case 'downarrow'
                    pow = 4^double(keyMod(hEvent, 'shift')); % 1 or 4
                    projScale = hFigProj.figData.boundScale*sqrt(2)^pow;
                    rescaleFigProj(hFigProj, projScale, obj.hCfg);

                case 'leftarrow' % go down one channel
                    if obj.projSites(1) > 1
                        obj.projSites = obj.projSites - 1;
                        obj.updateFigProj(false);
                    end

                case 'rightarrow' % go up one channel
                    if obj.projSites(end) < numel(obj.hCfg.siteMap)
                        obj.projSites = obj.projSites + 1;
                        obj.updateFigProj(false);
                    end

                case 'b' % background spikes
                    hFigProj.toggleVisible('background');

                case 'f' % toggle feature display
                    if strcmp(obj.hCfg.dispFeature, 'vpp')
                        obj.updateProjection(obj.hCfg.clusterFeature);
                    else
                        obj.updateProjection('vpp');
                    end

                case 'r' %reset view
                    obj.updateFigProj(true);

                case 's' %split
                    disp('Split: not implemented yet');
%                     if numel(obj.selected) == 1
%                         iCluster = obj.selected(1);
% 
%                         hFigProj.addPlot('hPoly', @impoly)
%                         polyPos = hFigProj.plotApply('hPoly', @getPosition);
% 
%                         XData = hFigProj.plotApply('foreground', @get, 'XData');
%                         YData = hFigProj.plotApply('foreground', @get, 'YData');
% 
%                         retained = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));
%                         jSites = unique(floor(XData(retained))) + 1;
%                         iSites = unique(floor(YData(retained))) + 1;
% 
%                         % return here
%                         hFigProj.addPlot('hSplit', @line, XData(retained), YData(retained), ...
%                                          'Color', obj.hCfg.mrColor_proj(3, :), ...
%                                          'Marker', '.', 'LineStyle', 'none');
% 
%                         dlgAns = questdlg('Split?', 'Confirmation', 'No');
% 
%                         hFigProj.rmPlot('hPoly');
%                         hFigProj.rmPlot('hSplit');
% 
%                         if strcmp(dlgAns, 'Yes')
%                             obj.splitCluster(iCluster, retained);
%                         end
%                     end

                case 'm' % merge clusters
                    hFigProj.wait(true);
                    obj.mergeSelected();
                    hFigProj.wait(false);

                case 'p' % toggle PCi v. PCj
                    if strcmp(obj.hCfg.dispFeature, 'pca')
                        % [1, 2] => [1, 3] => [2, 3] => [1, 2] => ...
                        obj.hCfg.pcPair = sort(mod(obj.hCfg.pcPair + 1, 3) + 1);
                        obj.updateFigProj(false);
                    end

                case 'h' % help
                    msgbox_(hFigProj.figData.helpText, 1);
            end % switch
        end

        function keyPressFigTime(obj, hObject, hEvent)
            %KEYPRESSFIGTIME Handle callbacks for keys pressed in time view
            hFigTime = obj.hFigs('FigTime');

            nSites = numel(obj.hCfg.siteMap);
            switch hEvent.Key
                case 'leftarrow' % go down one channel
                    factor = 3*double(keyMod(hEvent, 'shift')) + 1; % 1 or 4
                    obj.currentSite = max(obj.currentSite - factor, 1);
                    obj.updateFigTime(false);

                case 'rightarrow' % go up one channel
%                     if ~isVisible_(S_fig.hAx)
%                         msgbox_('Channel switching is disabled in the position view'); return;
%                     end
                    factor = 3*double(keyMod(hEvent, 'shift')) + 1; % 1 or 4
                    obj.currentSite = min(obj.currentSite + factor, nSites);
                    obj.updateFigTime(false);

                case 'uparrow'
                    pow = -(4^double(keyMod(hEvent, 'shift'))); % -1 or -4
                    rescaleFigTime(hFigTime, sqrt(2)^pow);
                    
                case 'downarrow' % change amp
                    pow = 4^double(keyMod(hEvent, 'shift')); % 1 or 4
                    rescaleFigTime(hFigTime, sqrt(2)^pow);

                case 'b' % toggle background spikes
                    hFigTime.figData.doPlotBG = hFigTime.toggleVisible('background');
                   
                case 'c' % compare pca across channels
                    disp('channel pca');
%                     hMsg = msgbox_('Plotting...');
%                     figure; hold on;
%                     [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(obj.selected(1));
%                     [~, mrPv1] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
%                     mrPv1 = norm_mr_(mrPv1);
% 
%                     if keyMod(event, 'control') %show chain of clusters
%                         trPv1 = mrPv1;
%                         iClu_next = get_next_clu_(S_clu, obj.selected(1));
%                         viClu_track = obj.selected(1);
%                         while ~isempty(iClu_next)
%                             [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(iClu_next);
%                             [~, mrPv1a] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
%                             mrPv1a = norm_mr_(mrPv1a);
%                             mrPv1 = flip_prinvec_(mrPv1a, mean(trPv1,3));
%                             trPv1 = cat(3, trPv1, mrPv1);
%                             viClu_track(end+1) = iClu_next;
% 
%                             iClu_next = get_next_clu_(S_clu, iClu_next);
%                         end
%                         multiplot(plot(nan,nan,'k'), 1, 1:size(trPv1,1), trPv1);
%             %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
%                         vcTitle = sprintf('PCA across chan: Clu %s', sprintf('%d,', viClu_track));
%                     elseif ~isempty(S0.iCluPaste)
%                         [mrWav_mean2, viSite1] = mrWav_int_mean_clu_(S0.iCluPaste);
%                         [~, mrPv2] = pca(mrWav_mean2, 'NumComponents', P.nPc_dip);
%                         mrPv2 = match_mrPv_(mrPv2, mrPv1);
%             %             mrPv2 = flip_prinvec_(mrPv2, mrPv1);
%                         mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
%                         mr2plot(norm_mr_(mrPv2), 'scale', 1, 'LineStyle', 'r--');
%                         vcTitle = sprintf('PCA across chan: Clu %d vs %d', obj.selected(1), S0.iCluPaste);
%                     else
%                         mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'r');
%                         vcTitle = sprintf('PCA across chan: Clu %d', obj.selected(1));
%                     end
%             %         mr2plot(mrPv1, 'scale', 1, 'LineStyle', 'k');
%                     grid on;
%                     title_(vcTitle);
%             %         if ~isempty(S0.iCluPaste)
%             %             compare_interp_(Sclu, obj.selected(1), S0.iCluPaste);
%             %         end
%                     try close(hMsg); catch; end

                case 'e' % export selected to workspace
                    disp('export');
%                     eval(sprintf('mrFet_clu%d = getDispFeaturesCluster(obj.selected(1));', obj.selected(1)));
%                     mrDist1 = squareform(pdist(mrFet1'));
%                     vrFet1 = sqrt(sum(mrFet1.^2));
%                     mrDist1 = bsxfun(@rdivide, mrDist1, vrFet1); %norm
%                     eval(sprintf('assignWorkspace_(mrFet_clu%d);', obj.selected(1)));

                case 'f' % toggle feature display
                    if strcmp(obj.hCfg.dispFeature, 'vpp')
                        obj.updateProjection(obj.hCfg.clusterFeature);
                    else
                        obj.updateProjection('vpp');
                    end

                case 'h' % help
                    msgbox_(hFigTime.figData.helpText, 1);

                case 'm' % merge
                    hFigTime.wait(true);
                    obj.mergeSelected();
                    hFigTime.wait(false);

                case 'r' % reset view
                    obj.updateFigTime(true);

                case 's' % split
                    if numel(obj.selected) == 1
                        iCluster = obj.selected(1);

                        hFigTime.addPlot('hPoly', @impoly)
                        polyPos = hFigTime.plotApply('hPoly', @getPosition);

                        XData = hFigTime.plotApply('foreground', @get, 'XData');
                        YData = hFigTime.plotApply('foreground', @get, 'YData');

                        retained = inpolygon(XData, YData, polyPos(:,1), polyPos(:,2));
                        hFigTime.addPlot('hSplit', @line, XData(retained), YData(retained), ...
                                         'Color', obj.hCfg.mrColor_proj(3, :), ...
                                         'Marker', '.', 'LineStyle', 'none');

                        dlgAns = questdlg('Split?', 'Confirmation', 'No');

                        hFigTime.rmPlot('hPoly');
                        hFigTime.rmPlot('hSplit');

                        if strcmp(dlgAns, 'Yes')
                            obj.splitCluster(iCluster, retained);
                        end
                    end
            end
        end

        function keyPressFigWav(obj, hObject, hEvent)
            %KEYPRESSFIGWAV Handle callbacks for keys pressed in main view
            hFigWav = obj.hFigs('FigWav');
            nSites = numel(obj.hCfg.siteMap);

            switch hEvent.Key
                case 'uparrow'
                    pow = -(4^double(keyMod(hEvent, 'shift'))); % -1 or -4
                    obj.maxAmp = rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, sqrt(2)^pow);
                    obj.updateCursorFigWav();

                case 'downarrow'
                    pow = 4^double(keyMod(hEvent, 'shift')); % 1 or 4
                    obj.maxAmp = rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, sqrt(2)^pow);
                    obj.updateCursorFigWav();

                case 'leftarrow' % select previous cluster
                    if keyMod(hEvent, 'shift')
                        selected_ = [obj.selected(1), max(obj.selected(end)-1, 1)];
                    else
                        selected_ = max(obj.selected(1)-1, 1);
                    end
                    obj.updateSelect(selected_);

                case 'rightarrow' % select next cluster
                    if keyMod(hEvent, 'shift')
                        selected_ = [obj.selected(1), min(obj.selected(end)+1, obj.hClust.nClusters)];
                    else
                        selected_ = min(obj.selected(1)+1, obj.hClust.nClusters);
                    end
                    obj.updateSelect(selected_);

                case 'home' % select first cluster
                    obj.updateSelect(1);
                    obj.keyPressFigWav([], struct('Key', 'z')); % zoom in

                case 'end' % select last cluster
                    obj.updateSelect(obj.hClust.nClusters);
                    obj.keyPressFigWav([], struct('Key', 'z')); % zoom in

                case {'d', 'backspace', 'delete'}
                    hFigWav.wait(true);
                    obj.deleteClusters();
                    hFigWav.wait(false);

                case 'm' % merge clusters
                    hFigWav.wait(true);
                    obj.mergeSelected();
                    hFigWav.wait(false);

                case 'n' % toggle spike count in clusters
                    obj.hCfg.fText = ~obj.hCfg.fText;
                    setFigWavXTicks(hFigWav, obj.hClust, obj.hCfg.fText);

                case 'space' % select most similar to currently selected
                    simScore = obj.hClust.simScore;
                    simScore(obj.selected(1), obj.selected(1)) = -inf;
                    [~, nextBest] = max(simScore(:, obj.selected(1)));
                    obj.updateSelect([obj.selected(1), nextBest]);

                case 's' % split
                    hFigWav.wait(true);
                    obj.autoSplit(true);
                    hFigWav.wait(false);

                case 'r' %reset view
                    hFigWav.wait(true);
                    hFigWav.axis([0, obj.hClust.nClusters + 1, 0, numel(obj.hCfg.siteMap) + 1]);
                    hFigWav.wait(false);

                case 'w' % toggle individual spike waveforms
                    hFigWav.toggleVisible('hSpkAll');

                case 'z' % zoom
                    if isempty(obj.selected)
                        obj.updateSelect(1);
                    else
                        iCluster = obj.selected(1);
                        iSite = obj.hClust.clusterSites(iCluster);
                        hFigWav.setWindow(iCluster + [-1, 1]*6, iSite + [-1, 1]*(obj.hCfg.maxSite*2+1), [0 obj.hClust.nClusters+1], [0 nSites+1]);
                    end

%                 case 'a', update_spikes_(S0); clu_info_(S0);
                case 'h'
                    msgbox_(hFigWav.figData.helpText, 1);

                case {'0', 'numpad0'}
                    obj.annotateUnit('to_delete', false); % TW

                case {'1', 'numpad1'}
                    obj.annotateUnit('single', false); % TW

                case {'2', 'numpad2'}
                    obj.annotateUnit('multi', false); % TW
%                 case 'p' %PSTH plot
%                     if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
%                     plot_raster_(S0, 1);
                otherwise
                    hFigWav.wait(false); %stop waiting
            end

            hFigWav.toForeground(); %change the focus back to the current object
        end

        function mouseClickFigSim(obj, xyPos, clickType)
            %MOUSECLICKFIGSIM Handle callbacks for mouse clicks in sim view
            xyPos = max(round(xyPos), [1 1]);
            if strcmp(clickType, 'normal') % left click
                obj.updateSelect(xyPos);
            end
        end

        function mouseClickFigWav(obj, xyPos, clickType)
            %MOUSECLICKFIGWAV Handle callbacks for mouse clicks in main view
            iCluster = round(xyPos(1)); % floor of x position
            if iCluster < 1 || iCluster > obj.hClust.nClusters
                return;
            end

            if strcmp(clickType, 'normal')
                obj.updateSelect(iCluster);
            elseif strcmp(clickType, 'alt') % right click, select secondary cluster
                obj.updateSelect([obj.selected(1) iCluster]);
            else                            % middle click, ignore
                disp(clickType);
                return;
            end

            hFig = obj.hFigs('FigWav');
            hFig.wait(true);
            % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z'
            % plot_raster_(S0);
            hFig.wait(false);
        end
    end

    %% SPECIFIC FIGURE METHODS
    methods
        function iCluster = plotSelectedWaveforms(obj, iCluster, plotKey)
            %PLOTSELECTEDWAVEFORMS Having selected a cluster, plot an
            %overlay on its waveforms
            hFigWav = obj.hFigs('FigWav');

            if strcmp(plotKey, 'selected2')
                colorMap = obj.hCfg.mrColor_proj(3, :); % red
            elseif strcmp(plotKey, 'selected1')
                colorMap = obj.hCfg.mrColor_proj(2, :); % black
            else
                return; % maybe some day we support selected3, 4, ...
            end

            if ~hFigWav.hasPlot(plotKey)
                hFigWav.addPlot(plotKey, nan, nan, 'Color', colorMap, 'LineWidth', 2);
            end

            if obj.hCfg.showRaw
                meanWf = obj.hClust.meanWfGlobalRaw(:, :, iCluster);
            else
                meanWf = obj.hClust.meanWfGlobal(:, :, iCluster);
            end

            % this is a hack until we hunt down and destroy all calls to
            % `multiplot` which are actually rescales (nargin <= 2)
%             userData = hFigWav.plotGet(plotKey, 'UserData');
%             if isempty(userData)
%                 userData = struct('maxAmp', obj.hCfg.maxAmp);
%             elseif ~isfield(userData, 'maxAmp')
%                 userData.maxAmp = obj.hCfg.maxAmp;
%             end 

            hFigWav.multiplot(plotKey, obj.maxAmp, getXRange(iCluster, obj.hCfg), meanWf);
            hFigWav.plotApply(plotKey, @uistack, 'top');
        end

        function updateCursorFigSim(obj)
            %UPDATECURSORFIGSIM
            if isempty(obj.selected) || ~obj.hasFig('FigSim')
                return;
            end

            iCluster = obj.selected(1);
            if numel(obj.selected) == 1
                jCluster = iCluster;
            else
                jCluster = obj.selected(2);
            end

            hFigSim = obj.hFigs('FigSim');

            % update crosshair cursor
            hFigSim.updatePlot('hCursorV', iCluster*[1, 1], 0.5 + [0, obj.hClust.nClusters]);
            if iCluster == jCluster
                colorH = obj.hCfg.mrColor_proj(2, :); % black
            else
                colorH = obj.hCfg.mrColor_proj(3, :); % red
            end
            hFigSim.updatePlot('hCursorH', 0.5 + [0, obj.hClust.nClusters], jCluster*[1, 1]);
            hFigSim.plotApply('hCursorH', @set, 'Color', colorH);

            % center on this pair of clusters
            hFigSim.axApply(@set, 'XLim', jrclust.utils.trimLim(iCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));
            hFigSim.axApply(@set, 'YLim', jrclust.utils.trimLim(jCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));

            scoreij = obj.hClust.simScore(iCluster, jCluster);
            hFigSim.axApply(@title, sprintf('Cluster %d vs. Cluster %d: %0.3f', iCluster, jCluster, scoreij), 'Interpreter', 'none', 'FontWeight', 'normal');
        end

        function updateCursorFigWav(obj)
            %UPDATECURSORFIGWAV Plot mean waveforms on top of selected
            %cluster(s)
            if isempty(obj.selected) || ~obj.hasFig('FigWav')
                return;
            end

            hFig = obj.hFigs('FigWav');
            if numel(obj.selected) == 1 % we've selected just one cluster (primary)
                hFig.rmPlot('selected2'); % if already selected, hide it
                obj.plotSelectedWaveforms(obj.selected, 'selected1');
            else % just plot #2 for now
                obj.plotSelectedWaveforms(obj.selected(1), 'selected1');
                obj.plotSelectedWaveforms(obj.selected(2), 'selected2');
            end
        end

        function updateFigCorr(obj)
            %UPDATEFIGCORR Plot cross correlation
            if isempty(obj.selected) || ~obj.hasFig('FigCorr')
                return;
            end

            doPlotFigCorr(obj.hFigs('FigCorr'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigHist(obj)
            %UPDATEFIGHIST Plot ISI histogram
            if isempty(obj.selected) || ~obj.hasFig('FigHist')
                return;
            end

            doPlotFigHist(obj.hFigs('FigHist'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigISI(obj)
            %UPDATEFIGISI Plot return map
            if isempty(obj.selected) || ~obj.hasFig('FigISI')
                return;
            end

            doPlotFigISI(obj.hFigs('FigISI'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigMap(obj)
            %UPDATEFIGMAP Plot probe map
            if isempty(obj.selected) || ~obj.hasFig('FigMap')
                return;
            end

            doPlotFigMap(obj.hFigs('FigMap'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigPos(obj)
            %UPDATEFIGPOS Plot cluster position on probe
            if isempty(obj.selected) || ~obj.hasFig('FigPos')
                return;
            end

            doPlotFigPos(obj.hFigs('FigPos'), obj.hClust, obj.hCfg, obj.selected, obj.maxAmp);
        end

        function updateFigProj(obj, doAutoscale)
            %UPDATEFIGPROJ
            if ~obj.hasFig('FigProj')
                return;
            end

            if nargin < 2
                doAutoscale = true;
            end

            hFigProj = obj.hFigs('FigProj');

            if isempty(hFigProj.figData) || ~isfield(hFigProj.figData, 'boundScale')
                boundScale = obj.maxAmp;
            else
                boundScale = hFigProj.figData.boundScale;
            end
            doPlotFigProj(hFigProj, obj.hClust, obj.projSites, obj.selected, boundScale);
            if doAutoscale
                autoScaleFigProj(hFigProj, obj.hClust, obj.selected);
            end
        end

        function updateFigRD(obj)
            %UPDATEFIGRD
            if ~obj.hasFig('FigRD')
                return;
            end

            doPlotFigRD(obj.hFigs('FigRD'), obj.hClust, obj.hCfg);
        end

        function updateFigSim(obj)
            %UPDATEFIGSIM
            if ~obj.hasFig('FigSim')
                return;
            end

            doPlotFigSim(obj.hFigs('FigSim'), obj.hClust, obj.hCfg);
        end

        function updateFigTime(obj, doAutoscale)
            %UPDATEFIGTIME
            if ~obj.hasFig('FigTime')
                return;
            end

            if nargin < 2
                doAutoscale = true;
            end

            hFigTime = obj.hFigs('FigTime');
            doPlotFigTime(hFigTime, obj.hClust, obj.hCfg, obj.selected, obj.maxAmp, obj.currentSite);
            if doAutoscale
                autoScaleFigTime(hFigTime, obj.hClust, obj.selected);
            end
        end

        function updateFigWav(obj)
            %UPDATEFIGWAV
            if ~obj.hasFig('FigWav')
                return;
            end

            doPlotFigWav(obj.hFigs('FigWav'), obj.hClust, obj.hCfg, obj.maxAmp);
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function addMenu(obj, hFig)
            drawnow;

            outerPosition = hFig.outerPosition;
            hFig.figSet('MenuBar','None');

            obj.hMenus('FileMenu') = hFig.uimenu('Label', 'File');
            uimenu(obj.hMenus('FileMenu'), 'Label', 'Save', 'Callback', @obj.saveFiles); % save_manual_
            uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .fig', 'Callback', @(hO, hE) obj.saveFigures('.fig')); % save_figures_
            uimenu(obj.hMenus('FileMenu'), 'Label', 'Save figures as .png', 'Callback', @(hO, hE) obj.saveFigures('.png')); % save_figures_
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export units to csv', 'Callback', @export_csv_, 'Separator', 'on');
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export unit qualities to csv', 'Callback', @(hO, hE)export_quality_);
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export selected mean unit waveforms', 'Callback', @(hO, hE)export_mrWav_clu_);
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export all waveforms from the selected unit', 'Callback', @(hO, hE)export_tnWav_spk_);
%             uimenu(obj.hMenus('FileMenu'), 'Label', 'Export firing rate for all units', 'Callback', @(hO, hE)export_rate_);
            uimenu(obj.hMenus('FileMenu'), 'Label', 'Exit', 'Callback', @(hO, hE) obj.endSession(), 'Separator', 'on', 'Accelerator', 'Q');

            obj.hMenus('EditMenu') = hFig.uimenu('Label', 'Edit');
            uimenu(obj.hMenus('EditMenu'), 'Label', '[M]erge', 'Callback', @(hO, hE) obj.mergeSelected());
%             uimenu(obj.hMenus('EditMenu'),'Label', 'Merge auto', 'Callback', @(hO, hE) merge_auto_());
            uimenu(obj.hMenus('EditMenu'), 'Label', '[D]elete', 'Callback', @(hO, hE) obj.deleteClusters(), 'Separator', 'on');
%             uimenu(obj.hMenus('EditMenu'),'Label', 'Delete auto', 'Callback', @(hO, hE) delete_auto_());
            uimenu(obj.hMenus('EditMenu'), 'Label', 'Delete annotated', 'Callback', @(hO, hE) obj.deleteAnnotated()); % TW
            uimenu(obj.hMenus('EditMenu'), 'Label', '[S]plit', 'Callback', @(hO, hE) obj.autoSplit(true), 'Separator', 'on');
            uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split max-chan', 'Callback', @(hO, hE) obj.autoSplit(false));
            uimenu(obj.hMenus('EditMenu'), 'Label', 'Auto split multi-chan', 'Callback', @(hO, hE) obj.autoSplit(true));

            obj.hMenus('ViewMenu') = uimenu(hFig, 'Label', 'View');
%             uimenu(obj.hMenus('ViewMenu'),'Label', 'Show traces', 'Callback', @(hO, hE)traces_());
            uimenu(obj.hMenus('ViewMenu'),'Label', 'View all [R]', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'r')));
            uimenu(obj.hMenus('ViewMenu'), 'Label', '[Z]oom selected', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'z')));
            uimenu(obj.hMenus('ViewMenu'), 'Label', '[W]aveform (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'w')));
            uimenu(obj.hMenus('ViewMenu'),'Label', '[N]umbers (toggle)', 'Callback', @(hO, hE) obj.keyPressFigWav([], struct('Key', 'n')));
            uimenu(obj.hMenus('ViewMenu'),'Label', 'Show raw waveform', 'Callback', @(hO, hE) obj.toggleRaw(hO))
%             %uimenu(obj.hMenus('ViewMenu'),'Label', 'Threshold by sites', 'Callback', @(hO, hE)keyPressFcn_thresh_(hFig, 'n'));
            uimenu(obj.hMenus('ViewMenu'),'Label', 'Reset window positions', 'Callback', @(hO, hE) obj.resetPositions());
% 
            obj.hMenus('ProjMenu') = uimenu(hFig,'Label', 'Projection');
            uimenu(obj.hMenus('ProjMenu'), 'Label', 'vpp', 'Callback', @(hO, hE) obj.updateProjection('vpp'));
            uimenu(obj.hMenus('ProjMenu'), 'Label', 'pca', 'Callback', @(hO, hE) obj.updateProjection('pca'));
            uimenu(obj.hMenus('ProjMenu'), 'Label', 'ppca', 'Callback', @(hO, hE) obj.updateProjection('ppca'));
            % uimenu(obj.hMenus('ProjMenu'), 'Label', 'cov', 'Callback', @(hO, hE) obj.updateProjection('cov'));
% 
%             mh_plot = uimenu(hFig,'Label','Plot');
%             uimenu(mh_plot, 'Label', 'All unit firing rate vs. aux. input', 'Callback', @(hO, hE)plot_aux_rate_);
%             uimenu(mh_plot, 'Label', 'Selected unit firing rate vs. aux. input', 'Callback', @(hO, hE)plot_aux_rate_(1));
% 
            obj.hMenus('InfoMenu') = uimenu(hFig, 'Label', '', 'Tag', 'InfoMenu');
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Annotate unit', 'Callback', @(hO, hE) obj.annotateUnit('', true));
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Single unit', 'Callback', @(hO, hE) obj.annotateUnit('single', false), 'Accelerator', '1');
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Multi unit', 'Callback', @(hO, hE) obj.annotateUnit('multi', false), 'Accelerator', '2');
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Noise', 'Callback', @(hO, hE) obj.annotateUnit('noise', false));
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Clear annotation', 'Callback', @(hO, hE) obj.annotateUnit('', false));
            uimenu(obj.hMenus('InfoMenu'), 'Label', 'Equal to', 'Callback', @(hO, hE) obj.annotateUnit('=', true));
% 
%             mh_help = uimenu(hFig,'Label','Help');
%             uimenu(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);

            drawnow;
            hFig.outerPosition = outerPosition;
        end

        function annotateUnit(obj, note, doConfirm)
            %ANNOTATEUNIT Add a note to a cluster
            iCluster = obj.selected(1);

            if nargin < 2
                note = obj.hClust.clusterNotes{iCluster};
            elseif isempty(note) % clear annotation
                obj.hClust.addNote(iCluster, '');
            end

            % set equal to another cluster?
            if ~isempty(note) && strcmp(note, '=') && numel(obj.selected) == 2
                note = sprintf('=%d', obj.selected(2));
            elseif ~isempty(note) && strcmp(note, '=')
                msgbox('Right-click another unit to set equal to selected unit');
                return;
            end

            if doConfirm
                newNote = inputdlg(sprintf('Cluster %d', iCluster), 'Annotation', 1, {note});
                if ~isempty(newNote)
                    obj.hClust.addNote(iCluster, newNote{1});
                end
            else
                obj.hClust.addNote(iCluster, note);
            end

            obj.updateMenu();
        end

        function autoSplit(obj, fMultisite)
            %AUTOSPLIT
            if numel(obj.selected) > 1
                return;
            end

            if obj.hClust.clusterCounts(obj.selected) < 2
                msgbox('At least two spikes required for splitting');
                return;
            end

%             hBox = msgbox('Splitting... (this closes automatically)');
            iCluster = obj.selected(1);
            iSite = obj.hClust.clusterSites(iCluster);

            if fMultisite
                spikeSites = obj.hCfg.siteNeighbors(1:end - obj.hCfg.nSitesExcl, iSite);
            else
                spikeSites = iSite;
            end

            iSpikes = obj.hClust.spikesByCluster{iCluster};

            sampledSpikes = jrclust.utils.getSampledWindows(obj.hClust, iSpikes, spikeSites);
            sampledSpikes = jrclust.utils.filtTouV(sampledSpikes, obj.hCfg);
            sampledSpikes = reshape(sampledSpikes, [], size(sampledSpikes, 3));

            % get vpp of cluster spikes on current site (in FigTime)
            localSpikes = jrclust.utils.getSampledWindows(obj.hClust, iSpikes, obj.currentSite);
            localSpikes = squeeze(jrclust.utils.filtTouV(localSpikes, obj.hCfg)); % TW calculate amplitudes on the fly
            localVpp = max(localSpikes) - min(localSpikes); % TW calculate amplitudes on the fly

            clusterTimes = obj.hClust.spikeTimes(iSpikes);
            [retained, splitFeatures] = doAutoSplit(sampledSpikes, [clusterTimes localVpp'], 2, obj.hCfg); %TW

            % create split figure
            hFigSplit = jrclust.views.Figure('FigSplit', [.5 0 .5 1], 'Split', false, false);
            hFigSplit.addSubplot('pcPlots', 2, 2);

            while true
                splitOff = ~retained;
                % plot PC2 vs. PC1
                hFigSplit.subplotApply('pcPlots', 1, @plot, ...
                                       splitFeatures(retained, 1), splitFeatures(retained, 2), 'b.', ...
                                       splitFeatures(splitOff, 1), splitFeatures(splitOff, 2), 'r.');
                hFigSplit.subplotApply('pcPlots', 1, @xlabel, 'PC 1');
                hFigSplit.subplotApply('pcPlots', 1, @ylabel, 'PC 2');

                % plot PC2 vs. PC3
                hFigSplit.subplotApply('pcPlots', 2, @plot, ...
                                       splitFeatures(retained, 3), splitFeatures(retained, 2), 'b.', ...
                                       splitFeatures(splitOff, 3), splitFeatures(splitOff, 2), 'r.');
                hFigSplit.subplotApply('pcPlots', 2, @xlabel, 'PC 3');
                hFigSplit.subplotApply('pcPlots', 2, @ylabel, 'PC 2');

                % plot PC3 vs. PC1
                hFigSplit.subplotApply('pcPlots', 3, @plot, ...
                                       splitFeatures(retained, 1), splitFeatures(retained, 3), 'b.', ...
                                       splitFeatures(splitOff, 1), splitFeatures(splitOff, 3), 'r.');
                hFigSplit.subplotApply('pcPlots', 3, @xlabel, 'PC 1');
                hFigSplit.subplotApply('pcPlots', 3, @ylabel, 'PC 3');

                % plot mean waveforms
                yMin = min(reshape(sampledSpikes, 1, []));
                yMax = max(reshape(sampledSpikes, 1, []));
                meanSamp = jrclust.utils.subsample(1:size(sampledSpikes, 1), 1000);

                hFigSplit.subplotApply('pcPlots', 4, @plot, ...
                                       mean(sampledSpikes(meanSamp, retained), 2), 'b');
                hFigSplit.subplotApply('pcPlots', 4, @hold, 'on');
                hFigSplit.subplotApply('pcPlots', 4, @plot, ...
                                       mean(sampledSpikes(meanSamp, splitOff), 2), 'r');
                hFigSplit.subplotApply('pcPlots', 4, @ylim, [yMin yMax]);
                hFigSplit.subplotApply('pcPlots', 4, @hold, 'off');
                hFigSplit.subplotApply('pcPlots', 4, @title, 'Mean spike waveforms');

                % ask if we want to split
                dlgAns = questdlg('Split?', 'Confirm split', ...
                                  'Yes', 'No', 'Manual', ... % options
                                  'No');                     % default

                if strcmp(dlgAns, 'Yes')
                    hFigSplit.close();
                    break;
                elseif strcmp(dlgAns, 'No')
                    hFigSplit.close();
                    return;
                else % Manual
                    dlgAns = questdlg('Select projection', '', ...
                                      'PC2 vs. PC1', 'PC2 vs. PC3', 'PC3 vs. PC1', ... % options
                                      'PC2 vs. PC1');                                  % default

                    if strcmp(dlgAns,'PC2 vs. PC1')
                        spIndex = 1;
                        pcX = 1;
                        pcY = 2;
                    elseif strcmp(dlgAns,'PC2 vs. PC3')
                        spIndex = 2;
                        pcX = 3;
                        pcY = 2;
                    else % PC3 vs. PC1
                        spIndex = 3;
                        pcX = 1;
                        pcY = 3;
                    end
                    % clear axes
                    hFigSplit.subplotApply('pcPlots', spIndex, @cla);

                    [featuresX, featuresY] = deal(splitFeatures(:, pcX), splitFeatures(:, pcY));
                    hFigSplit.subplotApply('pcPlots', spIndex, @plot, featuresX, featuresY, 'k.');

                    % user draws a polygon around features to keep
                    hPoly = hFigSplit.subplotApply('pcPlots', spIndex, @impoly);
                    polyPos = getPosition(hPoly);
                    retained = inpolygon(featuresX, featuresY, polyPos(:, 1), polyPos(:, 2));
                end
            end
            hFigSplit.close();

            obj.splitCluster(iCluster, retained);
        end

        function closeFigures(obj)
            %CLOSEFIGURES Close all open figures
            if obj.hasFig('FigWav') && obj.hFigs('FigWav').isReady
                hFigWav = obj.hFigs('FigWav');
                remove(obj.hFigs, 'FigWav'); % to prevent infinite recursion!
                hFigWav.close(); % calls killFigWav
            end

            try
                obj.figApply(@(hFig) hFig.close());
            catch ME
                warning(ME.identifier, 'Could not close figures: %s', ME.message);
            end

            obj.hFigs = containers.Map();
        end

        function deleteAnnotated(obj)
            %DELETEANNOTATED Delete clusters which are annotated to_delete
            deleteMe = find(strcmp(obj.hClust.clusterNotes, 'to_delete'));
            if ~isempty(deleteMe)
                obj.deleteClusters(deleteMe);
            end
        end

        function deleteClusters(obj, deleteMe)
            %DELETECLUSTERS Delete clusters either specified or selected
            if nargin < 2 && numel(obj.selected) > 1
                return;
            elseif nargin < 2
                deleteMe = obj.selected(1);
            end

            success = obj.hClust.deleteClusters(deleteMe);
            if success
                % save the new clustering
                deleted = strjoin(arrayfun(@num2str, deleteMe, 'UniformOutput', false), ', ');
                commitMsg = sprintf('%s;delete;%s', datestr(now, 31), deleted);
                obj.hClust.commit(commitMsg);

                % replot
                obj.updateFigWav();
                obj.updateFigRD(); % centers changed, need replotting
                obj.updateFigSim();
                obj.updateSelect(min(deleteMe));
            end
        end

        function res = figApply(obj, hFun, varargin)
            %FIGAPPLY Apply a function to all figures
            if nargout == 0 % "Too many output arguments"
                cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs), varargin{:});
            else
                res = cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs), varargin{:});
            end
        end

        function killFigWav(obj, hFig, ~)
            %KILLFIGWAV Destroy the main figure, close all other figures
            if ~obj.isEnding
                obj.endSession(); % we'll be back
            end

            delete(hFig);
        end

        function mergeSelected(obj)
            %MERGESELECTED Merge a pair of clusters
            if numel(obj.selected) < 2
                return;
            end
            iCluster = min(obj.selected);
            jCluster = max(obj.selected);

            success = obj.hClust.mergeClusterPair(iCluster, jCluster);
            if success
                % save the new clustering
                commitMsg = sprintf('%s;merge;%d;%d', datestr(now, 31), iCluster, jCluster);
                obj.hClust.commit(commitMsg);

                % replot
                obj.updateFigWav();
                obj.updateFigRD(); % centers changed, need replotting
                obj.updateFigSim();
                obj.updateSelect(min(obj.selected));
            end
        end

        function plotAllFigures(obj)
            %PLOTALLFIGURES Plot all figures
            if isempty(obj.hFigs)
                obj.spawnFigures();
            elseif ~all(obj.figApply(@(hFig) hFig.isReady)) % clean up from an aborted session
                obj.closeFigures();
                obj.spawnFigures();
            end

            % plot sim score figure
            if obj.hasFig('FigSim')
                % set key and mouse handles
                hFigSim = obj.hFigs('FigSim');
                hFigSim.hFunKey = @obj.keyPressFigSim;
                hFigSim.setMouseable(@obj.mouseClickFigSim);
                obj.updateFigSim();
            end

%                 case 'p' %PSTH plot
%                 if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
%                 plot_raster_(S0, 1);

            % plot feature projection
            if obj.hasFig('FigProj')
                % set key and mouse handles
                hFigProj = obj.hFigs('FigProj');
                hFigProj.hFunKey = @obj.keyPressFigProj;
                hFigProj.setMouseable(); % no special mouse function
            end

            % plot amplitude vs. time
            if obj.hasFig('FigTime')
                % set key and mouse handles
                hFigTime = obj.hFigs('FigTime');
                hFigTime.hFunKey = @obj.keyPressFigTime;
                hFigTime.setMouseable(); % no special mouse function
            end

            % plot main waveform view
            if obj.hasFig('FigWav')
                % set key and mouse handles
                hFigWav = doPlotFigWav(obj.hFigs('FigWav'), obj.hClust, obj.hCfg, obj.maxAmp);

                hFigWav.hFunKey = @obj.keyPressFigWav;
                hFigWav.setMouseable(@obj.mouseClickFigWav);

                % make this guy the key log
                hFigWav.figSet('CloseRequestFcn', @obj.killFigWav);
                obj.addMenu(hFigWav);
            end

            % plot rho-delta figure
            obj.updateFigRD();

            % select first cluster (also plots other figures)
            obj.updateSelect(1);

            % zoom in on first cluster
            obj.keyPressFigWav([], struct('Key', 'z')); % zoom in
        end

        function spawnFigures(obj)
            %SPAWNFIGURES Create new figures
            obj.hFigs = doSpawnFigures(obj.hCfg);
        end

        function splitCluster(obj, iCluster, retained)
            %SPLITCLUSTER Split off a cluster given retained spikes
            iSpikes = obj.hClust.spikesByCluster{iCluster};

            [success, retained] = obj.hClust.splitCluster(iCluster, iSpikes(retained));
            if success
                % save the new clustering
                retained = strjoin(arrayfun(@num2str, retained, 'UniformOutput', false), ',');
                commitMsg = sprintf('%s;split;%d;%s', datestr(now, 31), iCluster, retained);
                obj.hClust.commit(commitMsg);

                % replot
                obj.updateFigWav();
                obj.updateFigSim();
                obj.updateSelect([iCluster, iCluster + 1]);
            end
        end

        function toggleRaw(obj, hMenu)
            %TOGGLERAW Toggle raw waveform display
            showRaw = ~obj.hCfg.showRaw;
            if obj.hasFig('FigWav')
                hFigWav = obj.hFigs('FigWav');
                hFigWav.wait(true);
            end

            if showRaw && isempty(obj.hClust.meanWfGlobalRaw)
                obj.hClust.computeMeanWaveforms([], true);
            end

            set(hMenu, 'Checked', jrclust.utils.ifEq(showRaw, 'on', 'off'));
            obj.hCfg.showRaw = showRaw;

            % replot
            obj.updateFigWav();
            obj.updateFigSim();
            obj.updateSelect(obj.selected);

            if obj.hasFig('FigWav')
                hFigWav = obj.hFigs('FigWav');
                hFigWav.wait(false);
            end
        end

        function updateMenu(obj)
            % update menu entry to indicate selected clusters
            if numel(obj.selected) > 1 && obj.hasMenu('InfoMenu')
                menuLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', obj.selected(1), ...
                                    obj.hClust.clusterNotes{obj.selected(1)}, obj.selected(2), ...
                                    obj.hClust.clusterNotes{obj.selected(2)});
                set(obj.hMenus('InfoMenu'), 'Label', menuLabel);
            elseif obj.hasMenu('InfoMenu')
                menuLabel = sprintf('Unit %d "%s"', obj.selected(1), obj.hClust.clusterNotes{obj.selected(1)});
                set(obj.hMenus('InfoMenu'), 'Label', menuLabel);
            end
        end

        function updateProjection(obj, proj)
            %UPDATEPROJECTION Update the feature projection
            try
                obj.hCfg.dispFeature = proj;
            catch
                obj.hCfg.dispFeature = 'vpp';
            end

            % set menu items checked or unchecked where appropriate
            hProjMenu = obj.hMenus('ProjMenu');
            for i = 1:numel(hProjMenu.Children)
                if strcmp(hProjMenu.Children(i).Text, obj.hCfg.dispFeature)
                    set(hProjMenu.Children(i), 'Checked', 'on');
                else
                    set(hProjMenu.Children(i), 'Checked', 'off');
                end
            end

            obj.updateFigProj(true);
            obj.updateFigTime(true);
        end

        function updateSelect(obj, iClusters)
            %UPDATESELECT Select a (pair of) cluster(s) across all views
            iClusters = min(max(iClusters, 1), obj.hClust.nClusters);
            if numel(iClusters) > 2
                iClusters = iClusters(1:2);
            elseif isempty(iClusters)
                iClusters = 1;
            end
            if numel(iClusters) == 2 && diff(iClusters) == 0
                iClusters = iClusters(1);
            end

            obj.selected = iClusters;

            iSite = obj.hClust.clusterSites(obj.selected(1));

            % update current site for amplitude view
            obj.currentSite = iSite;
            % limit the number of sites to display in the feature projection view
            nSites = min(obj.hCfg.nSitesFigProj, size(obj.hCfg.siteNeighbors, 1)); % by request

            % center sites around cluster center site
            if nSites < size(obj.hCfg.siteNeighbors, 1)
                obj.projSites = iSite:iSite + nSites - 1;
                if obj.projSites(end) > max(obj.hCfg.siteMap) % correct for overshooting
                    obj.projSites = obj.projSites - max(obj.projSites) + max(obj.hCfg.siteMap);
                end
            else
                obj.projSites = sort(obj.hCfg.siteNeighbors(:, iSite), 'ascend');
            end

            % update plots
            obj.updateFigCorr();
            obj.updateFigHist();
            obj.updateFigISI();
            obj.updateFigMap();
            obj.updateFigPos();
            obj.updateFigProj(true);
            obj.updateFigTime(true);

            % update cursors
            obj.updateCursorFigWav();
            obj.updateCursorFigSim();

            % update menu
            obj.updateMenu();

            %obj.keyPressFigWav([], struct('Key', 'z'));
        end
    end

    %% USER METHODS
    methods
        function allToForeground(obj)
            %ALLTOFOREGROUND Bring all figures to the foreground
            obj.figApply(@(hFig) hFig.toForeground());
        end

        function beginSession(obj)
            %BEGINSESSION Start curating clusters
            if ~isempty(obj.hFigs) % session already running
                return;
            end
            obj.isEnding = false;
            obj.cRes = struct('hClust', obj.hClust);
            obj.selected = 1;
            obj.currentSite = obj.hClust.clusterSites(1);
            obj.maxAmp = obj.hCfg.maxAmp;
            obj.plotAllFigures();
        end

        function res = endSession(obj)
            %ENDSESSION Finish curating and return results
            if nargout == 0
                obj.saveFiles(); % ask
            end
            obj.isEnding = true;
            obj.closeFigures();

            obj.currentSite = [];
            obj.selected = [];
            res = obj.cRes;
            obj.cRes = [];
            obj.isEnding = false; % ended
        end

        function hf = hasFig(obj, figKey)
            %HASFIG Return true if we have a figure by key
            hf = ischar(figKey) && isKey(obj.hFigs, figKey);
        end

        function hm = hasMenu(obj, menuKey)
            %HASMENU Return true if we have a menu item by key
            hm = ischar(menuKey) && isKey(obj.hMenus, menuKey);
        end

        function resetPositions(obj)
            obj.figApply(@(hFig) hFig.resetPos());
        end

        function saveFigures(obj, ext)
            fprintf('saving figs as %s\n', ext);
        end

        function saveFiles(obj)
            %SAVEFILES Save files to disk
            disp('saveFiles');
        end
    end

    %% GETTERS/SETTERS
    methods
        % hCfg
        function hc = get.hCfg(obj)
            if ~isempty(obj.hClust)
                hc = obj.hClust.hCfg;
            else
                hc = [];
            end
        end

        % hClust
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = obj.hClust;
            end
        end
    end
end

