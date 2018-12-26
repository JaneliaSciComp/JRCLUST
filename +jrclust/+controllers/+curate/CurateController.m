classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters

    properties (Access=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hCfg;           % Config object
        hClust;         % Clustering object
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
        function obj = CurateController(hCfg)
            obj.hCfg = hCfg;
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
        function keyPressFigSim(obj, hObject, hEvent)
            %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
            hFigSim = obj.hFigs('hFigSim');
            switch hEvent.Key
                case {'d', 'backspace', 'delete'} % delete
                    % ui_delete_(S0);
                    disp('delete');

                case 'k' % kilosort template similarity view
                    if obj.hCfg.fImportKsort
                        % plot_FigWavCor_(S0, 'simscore');
                        disp('kilosort');
                    end

                case 'm' % merge
                    % ui_merge_(S0);
                    disp('merge');

                case 's' % split
                    % auto_split_(1);
                    disp('split');

                case 'w' % waveform-based sim view
                    % if get_set_(S0.P, 'fImportKsort', 0)
                    %     plot_FigWavCor_(S0, 'wavecor');
                    % end
                    if obj.hCfg.fImportKsort
                        % plot_FigWavCor_(S0, 'simscore');
                        disp('waveform');
                    end
            end % switch
        end

        function keyPressFigProj(obj, hObject, hEvent)
            hFigProj = obj.hFigs('hFigProj');
            %S_plot1 = get(S_fig.hPlot1, 'UserData');
            %viSites_show = S_plot1.viSites_show;
            hFigProj.wait(true);

            switch lower(hEvent.Key)
                case {'uparrow', 'downarrow'}
                    %rescaleFigProj(hEvent, hObject, S_fig, S0);
                    disp('rescale');

                case {'leftarrow', 'rightarrow'} % change channels
                    disp('change channels');
%                     fPlot = 0;
%                     if strcmpi(hEvent.Key, 'leftarrow')
%                         if min(S_fig.viSites_show)>1
%                             S_fig.viSites_show=S_fig.viSites_show-1;
%                             fPlot = 1;
%                         end
%                     else
%                         if max(S_fig.viSites_show) < max(hCfg.viSite2Chan)
%                             S_fig.viSites_show=S_fig.viSites_show+1;
%                             fPlot = 1;
%                         end
%                     end
%                     if fPlot
%                         set(hObject, 'UserData', S_fig);
%                         S0.P.viSites_show = S_fig.viSites_show;
%                         plot_FigProj_(S0);
%                     end

                case 'r' %reset view
                    disp('reset view');
%                     axis_([0 numel(viSites_show) 0 numel(viSites_show)]);

                case 's' %split
                    disp('split');
%                     figure_wait_(0);
%                     if ~isempty(S0.iCluPaste)
%                         msgbox_('Select one cluster to split'); return;
%                     end
%                     S_plot1 = select_polygon_(S_fig.hPlot1);
%                     if ~isempty(S_plot1)
%                         [fSplit, vlIn] = plot_split_(S_plot1);
%                         if fSplit
%                             hClust = split_clu_(S0.iCluCopy, vlIn);
%                         else
%                             update_plot2_proj_();
%                             %                 delete_multi_(S_plot1.hPoly);
%                         end
%                     end

                case 'm'
                    disp('merge');
                    %ui_merge_(S0);

                case 'f'
                    disp('toggle features');
%                     if get_set_(hCfg, 'fImportKsort', 0) && strcmpi(hCfg.vcFet_show, 'vpp')
%                         hCfg.vcFet_show = 'kilosort';
%                     elseif strcmpi(hCfg.vcFet_show, 'vpp')
%                         hCfg.vcFet_show = hCfg.vcFet;
%                     else
%                         hCfg.vcFet_show = 'vpp';
%                     end
% 
%                     set0_(hCfg);
%                     try
%                         plot_FigProj_();
%                     catch
%                         hCfg.vcFet_show = 'vpp';
%                         set0_(hCfg);
%                         plot_FigProj_();
%                     end

                case 'p' % toggle PCi v. PCj
                    disp('toggle pc pair');
%                     if get_set_(hCfg, 'fImportKsort', 0) && strcmpi(hCfg.vcFet_show, 'kilosort')
%                         pcPair = get_set_(S0, 'pcPair', [1 2]);
% 
%                         if all(pcPair == [1 2])
%                             S0.pcPair = [1 3];
%                         elseif all(pcPair == [1 3])
%                             S0.pcPair = [2 3];
%                         else
%                             S0.pcPair = [1 2];
%                         end
% 
%                         set(0, 'UserData', S0);
% 
%                         plot_FigProj_(S0);
%                     end

                case 'b' % background spikes
                    disp('toggle background');
                    %toggleVisible_(S_fig.hPlot0);

                case 'h' % help
                    disp('halp');
                    %msgbox_(S_fig.csHelp, 1);
            end % switch

            hFigProj.wait(false);
        end

        function keyPressFigTime(obj, hObject, hEvent)
            %KEYPRESSFIGTIME Handle callbacks for keys pressed in time view
            hFigTime = obj.hFigs('hFigTime');
            nSites = numel(obj.hCfg.siteMap);

            switch hEvent.Key
                case {'leftarrow', 'rightarrow'} % switch channel
                    disp('change channel');
%                     if ~isVisible_(S_fig.hAx)
%                         msgbox_('Channel switching is disabled in the position view'); return;
%                     end
                    factor = key_modifier_(hEvent, 'shift')*3 + 1;
                    if strcmpi(hEvent.Key, 'rightarrow')
                        obj.currentSite = min(obj.currentSite + factor, nSites);
                    else
                        obj.currentSite = max(obj.currentSite - factor, 1);
                    end
                    obj.updateFigTime();

                case {'uparrow', 'downarrow'} % change amp
                    disp('change amp');
%                     if ~isVisible_(S_fig.hAx)
%                         msgbox_('Zoom is disabled in the position view'); return;
%                     end
%                     rescaleFigTime(event, S0, hCfg);

                case 'r' %reset view
                    disp('reset');
%                     if ~isVisible_(S_fig.hAx)
%                         return;
%                     end
%                     axis_(S_fig.hAx, [S_fig.time_lim, S_fig.vpp_lim]);
%                     imrect_set_(S_fig.hRect, S_fig.time_lim, S_fig.vpp_lim);

                case 'm' % merge
                    disp('merge');
%                     ui_merge_(S0);

                case 'h' % help
                    msgbox_(hFigTime.figData.csHelp, 1);

                case 'b' % toggle background spikes
                    disp('toggle background spikes');
%                     if isVisible_(S_fig.hAx)
%                         S_fig.doPlotBG = toggleVisible_(S_fig.hPlot0);
%                     else
%                         S_fig.doPlotBG = toggleVisible_(S_fig.hPlot0_track);
%                     end
%                     set(hFig, 'UserData', S_fig);

%                 case 'z' % track depth
%                     disp('FigTime:''z'' not implemented yet');
%                     plot_SpikePos_(S0, event);

                case 's' % split
                    disp('split');
%                     if ~isempty(S0.iCluPaste)
%                         msgbox_('Select one cluster'); return;
%                     end
%                     try
%                         hPoly = impoly_();
%                         if isempty(hPoly); return ;end
%                         mrPolyPos = getPosition(hPoly);
%                         vrX1 = double(get(S_fig.hPlot1, 'XData'));
%                         vrY1 = double(get(S_fig.hPlot1, 'YData'));
%                         vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
%                         hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
%                         if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'yes')
%                             split_clu_(S0.iCluCopy, vlIn);
%                         end
%                         delete_multi_(hPoly, hSplit);
%                     catch
%                         disp(lasterror());
%                     end

                case 'p' % update projection
                    disp('update projection');
%                     vrPos = getPosition(S_fig.hRect);
%                     tlim_proj = [vrPos(1), sum(vrPos([1,3]))];
%                     hCfg.tlim_proj = tlim_proj;
%                     plot_FigProj_(S0);

%                 case 'f' % feature display instead of amplitude display
%                     if strcmpi(P.vcFet_show, 'fet')
%                         P.vcFet_show = 'vpp';
%                     else
%                         P.vcFet_show = 'fet';
%                     end
%                     set0_(P);
%                     update_FigTime_();

                case 'c' % compare pca across channels
                    disp('channel pca');
%                     hMsg = msgbox_('Plotting...');
%                     figure; hold on;
%                     [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(S0.iCluCopy);
%                     [~, mrPv1] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
%                     mrPv1 = norm_mr_(mrPv1);
% 
%                     if key_modifier_(event, 'control') %show chain of clusters
%                         trPv1 = mrPv1;
%                         iClu_next = get_next_clu_(S_clu, S0.iCluCopy);
%                         viClu_track = S0.iCluCopy;
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
%                         vcTitle = sprintf('PCA across chan: Clu %d vs %d', S0.iCluCopy, S0.iCluPaste);
%                     else
%                         mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'r');
%                         vcTitle = sprintf('PCA across chan: Clu %d', S0.iCluCopy);
%                     end
%             %         mr2plot(mrPv1, 'scale', 1, 'LineStyle', 'k');
%                     grid on;
%                     title_(vcTitle);
%             %         if ~isempty(S0.iCluPaste)
%             %             compare_interp_(Sclu, S0.iCluCopy, S0.iCluPaste);
%             %         end
%                     try close(hMsg); catch; end
% 
%                 case 'f' %feature export
%                     eval(sprintf('mrFet_clu%d = getDispFeaturesCluster(S0.iCluCopy);', S0.iCluCopy));
%                     mrDist1 = squareform(pdist(mrFet1'));
%                     vrFet1 = sqrt(sum(mrFet1.^2));
%                     mrDist1 = bsxfun(@rdivide, mrDist1, vrFet1); %norm
%                     eval(sprintf('assignWorkspace_(mrFet_clu%d);', S0.iCluCopy));

                case 'e' % export selected to workspace
                    disp('export');
            end
        end

        function keyPressFigWav(obj, hObject, hEvent)
            %KEYPRESSFIGWAV Handle callbacks for keys pressed in main view
            hFigWav = obj.hFigs('hFigWav');
            nSites = numel(obj.hCfg.siteMap);

            switch hEvent.Key
                case 'uparrow'
                    pow = 4^double(any(strcmpi(hEvent.Modifier, 'shift')));
                    obj.maxAmp = rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, 1/sqrt(2)^pow);

                    obj.updateCursorFigWav();

                case 'downarrow'
                    pow = 4^double(any(strcmpi(hEvent.Modifier, 'shift')));
                    obj.maxAmp = rescaleFigWav(hFigWav, obj.hClust, obj.hCfg, obj.maxAmp, sqrt(2)^pow);

                    obj.updateCursorFigWav();

                    %rescale_FigWav_(hEvent, S0, hCfg);
                    %clu_info_(S0); %update figpos

                case 'leftarrow' % select previous cluster
                    if any(strcmpi(hEvent.Modifier, 'shift'))
                        selected_ = [obj.selected(1), max(obj.selected(end)-1, 1)];
                    else
                        selected_ = max(obj.selected(1)-1, 1);
                    end
                    obj.updateSelect(selected_);

                case 'rightarrow' % select next cluster
                    if any(strcmpi(hEvent.Modifier, 'shift'))
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
%                     S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0); %select first clu
%                     if strcmpi(hEvent.Key, 'home') || strcmpi(hEvent.Key, 'end') %'z' to recenter
%                         S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0);
%                     end
                case 'm' % merge clusters
                    obj.mergeSelected();

                case 'space' % select most similar to currently selected
                    simScore = obj.hClust.simScore;
                    simScore(obj.selected(1), obj.selected(1)) = -inf;
                    [~, nextBest] = max(simScore(:, obj.selected(1)));
                    obj.updateSelect([obj.selected(1), nextBest]);

                case 's' % split
                    % auto_split_(1, S0);
                    disp('split');

                case 'r' %reset view
                    %figure_wait_(1);
                    %axis_([0, S0.S_clu.nClusters + 1, 0, numel(hCfg.viSite2Chan) + 1]);
                    %figure_wait_(0);
                    disp('reset');

                case {'d', 'backspace', 'delete'}
                    obj.deleteClusters();
                    disp('delete');

                case 'z' % zoom
                    if isempty(obj.selected)
                        obj.updateSelect(1);
                    else
                        iCluster = obj.selected(1);
                        iSite = obj.hClust.clusterSites(iCluster);
                        hFigWav.setWindow(iCluster + [-1, 1]*6, iSite + [-1, 1]*(obj.hCfg.maxSite*2+1), [0 obj.hClust.nClusters+1], [0 nSites+1]);
                    end

%                 case 'c', plot_FigCorr_(S0);
%                 case 'v', plot_FigIsi_(S0);
%                 case 'a', update_spikes_(S0); clu_info_(S0);
%                 case 'f', clu_info_(S0);
%                 case 'h', msgbox_(figData.csHelp, 1);
%                 case '0', unit_annotate_([],[], 'to_delete'); % TW
%                 case '1', unit_annotate_([],[], 'single'); % TW
%                 case '2', unit_annotate_([],[], 'multi'); % TW
%                 case 'numpad0', unit_annotate_([],[], 'to_delete'); % TW
%                 case 'numpad1', unit_annotate_([],[], 'single'); % TW
%                 case 'numpad2', unit_annotate_([],[], 'multi'); % TW
                case 'w'
                    hFigWav.toggleVisible('hSpkAll'); % toggle spike waveforms
%                 case 'n'
%                 fText = get_set_(figData, 'fText', get_set_(hCfg, 'fText', 1));
%                 setFigWavXTicks(figData, hClust, ~fText);
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

                % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
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

            hFig = obj.hFigs('hFigWav');
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
            hFigWav = obj.hFigs('hFigWav');

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
            hFigWav.uistack(plotKey, 'top');
        end

        function updateCursorFigSim(obj)
            %UPDATECURSORFIGSIM
            if isempty(obj.selected) || ~obj.hasFig('hFigSim')
                return;
            end

            iCluster = obj.selected(1);
            if numel(obj.selected) == 1
                jCluster = iCluster;
            else
                jCluster = obj.selected(2);
            end

            hFigSim = obj.hFigs('hFigSim');

            % update crosshair cursor
            hFigSim.updatePlot('hCursorV', iCluster*[1, 1], 0.5 + [0, obj.hClust.nClusters]);
            if iCluster == jCluster
                colorH = obj.hCfg.mrColor_proj(2, :); % black
            else
                colorH = obj.hCfg.mrColor_proj(3, :); % red
            end
            hFigSim.updatePlot('hCursorH', 0.5 + [0, obj.hClust.nClusters], jCluster*[1, 1]);
            hFigSim.plotSet('hCursorH', 'Color', colorH);

            % center on this pair of clusters
            hFigSim.axSet('XLim', jrclust.utils.trimLim(iCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));
            hFigSim.axSet('YLim', jrclust.utils.trimLim(jCluster + [-6, 6], 0.5 + [0, obj.hClust.nClusters]));

            scoreij = obj.hClust.simScore(iCluster, jCluster);
            hFigSim.title(sprintf('Cluster %d vs. Cluster %d: %0.3f', iCluster, jCluster, scoreij));
        end

        function updateCursorFigWav(obj)
            %UPDATECURSORFIGWAV
            if isempty(obj.selected) || ~obj.hasFig('hFigWav')
                return;
            end

            hFig = obj.hFigs('hFigWav');
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
            if isempty(obj.selected) || ~obj.hasFig('hFigCorr')
                return;
            end

            % hFigCorr = doPlotFigCorr(obj.hFigs('hFigCorr'), obj.hClust, obj.hCfg, obj.selected);
            % obj.hFigs('hFigCorr') = hFigCorr;
            doPlotFigCorr(obj.hFigs('hFigCorr'), obj.hClust, obj.hCfg, obj.selected);
        end
        %                 case 'i', plot_FigHist_(S0); %ISI histogram

        function updateFigHist(obj)
            %UPDATEFIGHIST Plot ISI histogram
            if isempty(obj.selected) || ~obj.hasFig('hFigHist')
                return;
            end
            
            % hFigHist = doPlotFigHist(obj.hFigs('hFigHist'), obj.hClust, obj.hCfg, obj.selected);
            % obj.hFigs('hFigHist') = hFigHist;
            doPlotFigHist(obj.hFigs('hFigHist'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigISI(obj)
            %UPDATEFIGISI Plot return map
            if isempty(obj.selected) || ~obj.hasFig('hFigISI')
                return;
            end

            % hFigISI = doPlotFigISI(obj.hFigs('hFigISI'), obj.hClust, obj.hCfg, obj.selected);
            % obj.hFigs('hFigISI') = hFigISI;
            doPlotFigISI(obj.hFigs('hFigISI'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigMap(obj)
            %UPDATEFIGMAP Plot probe map
            if isempty(obj.selected) || ~obj.hasFig('hFigMap')
                return;
            end

            % hFigMap = doPlotFigMap(obj.hFigs('hFigMap'), obj.hClust, obj.hCfg, obj.selected);
            % obj.hFigs('hFigMap') = hFigMap;
            doPlotFigMap(obj.hFigs('hFigMap'), obj.hClust, obj.hCfg, obj.selected);
        end

        function updateFigPos(obj)
            %UPDATEFIGPOS Plot cluster position on probe
            if isempty(obj.selected) || ~obj.hasFig('hFigPos')
                return;
            end

            % hFigPos = doPlotFigPos(obj.hFigs('hFigPos'), obj.hClust, obj.hCfg, obj.selected, obj.maxAmp);
            % obj.hFigs('hFigPos') = hFigPos;
            doPlotFigPos(obj.hFigs('hFigPos'), obj.hClust, obj.hCfg, obj.selected, obj.maxAmp);
        end

        function updateFigProj(obj)
            %UPDATEFIGPROJ
            iSite = obj.hClust.clusterSites(obj.selected(1));

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

            doPlotFigProj(obj.hFigs('hFigProj'), obj.hClust, obj.projSites, obj.selected, obj.maxAmp);
        end

        function updateFigTime(obj)
            %UPDATEFIGTIME
            if ~obj.hasFig('hFigTime')
                return;
            end

            hFigTime = obj.hFigs('hFigTime');

            % plot background spikes
            [bgFeatures, bgTimes, yLabel] = getFigTimeFeatures(obj.hClust, obj.currentSite);
            hFigTime.toggleVisible('background', hFigTime.figData.doPlotBG);
            hFigTime.updatePlot('background', bgTimes, bgFeatures);

            % plot foreground spikes
            [fgFeatures, fgTimes] = getFigTimeFeatures(obj.hClust, obj.currentSite, obj.selected(1));
            hFigTime.updatePlot('foreground', fgTimes, fgFeatures);

            % plot secondary foreground spikes
            if numel(obj.selected) == 2
                [fgFeatures2, fgTimes2] = getFigTimeFeatures(obj.hClust, obj.currentSite, obj.selected(2));
                hFigTime.updatePlot('foreground2', fgTimes2, fgFeatures2);
            else % or hide them
                hFigTime.hidePlot('foreground2');
            end

            hFigTime.axSet('YLim', [0, 1] * obj.maxAmp);
            imrect_set_(hFigTime, 'hRect', [], [0, 1] * obj.maxAmp);
            hFigTime.grid('on');
            hFigTime.ylabel(yLabel);
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function addMenu(obj, hFig)
            drawnow;

            outerPosition = hFig.outerPosition;
            hFig.figSet('MenuBar','None');

            obj.hMenus('hMenuFile') = hFig.uimenu('Label', 'File');
            uimenu(obj.hMenus('hMenuFile'), 'Label', 'Save', 'Callback', @obj.saveFiles); % save_manual_
            uimenu(obj.hMenus('hMenuFile'), 'Label', 'Save figures as .fig', 'Callback', @(h, e) obj.saveFigures('.fig')); % save_figures_
            uimenu(obj.hMenus('hMenuFile'), 'Label', 'Save figures as .png', 'Callback', @(h, e) obj.saveFigures('.png')); % save_figures_
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Describe', 'Callback', @(h,e) msgbox_(describe_()), 'Separator', 'on');
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Edit prm file', 'Callback', @edit_prm_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Reload prm file', 'Callback', @reload_prm_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export units to csv', 'Callback', @export_csv_, 'Separator', 'on');
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export unit qualities to csv', 'Callback', @(h,e)export_quality_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export selected mean unit waveforms', 'Callback', @(h,e)export_mrWav_clu_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export all waveforms from the selected unit', 'Callback', @(h,e)export_tnWav_spk_);
%             uimenu(obj.hMenus('hMenuFile'), 'Label', 'Export firing rate for all units', 'Callback', @(h,e)export_rate_);
            uimenu(obj.hMenus('hMenuFile'), 'Label', 'Exit', 'Callback', @(h, e) obj.endSession, 'Separator', 'on', 'Accelerator', 'Q');

%             mh_edit = uimenu(hFig,'Label','Edit');
%             uimenu(mh_edit,'Label', '[M]erge', 'Callback', @(h,e) keyPressFcn_cell_(hFig, 'm'));
%             uimenu(mh_edit,'Label', 'Merge auto', 'Callback', @(h,e) merge_auto_());
%             uimenu(mh_edit,'Label', '[D]elete', 'Callback', @(h,e) keyPressFcn_cell_(hFig, 'd'), 'Separator', 'on');
%             uimenu(mh_edit,'Label', 'Delete auto', 'Callback', @(h,e) delete_auto_());
%             uimenu(mh_edit,'Label', 'Delete annotated', 'Callback', @(h,e) delete_annotate()); % TW
%             uimenu(mh_edit,'Label', '[S]plit', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'), 'Separator', 'on');
%             uimenu(mh_edit,'Label', 'Auto split max-chan', 'Callback', @(h,e)auto_split_(0));
%             uimenu(mh_edit,'Label', 'Auto split multi-chan', 'Callback', @(h,e)auto_split_(1));
%             uimenu(mh_edit,'Label', 'Annotate', 'Callback', @(h,e)unit_annotate_());
% 
%             mh_view = uimenu(hFig, 'Label','View');
%             uimenu(mh_view,'Label', 'Show traces', 'Callback', @(h,e)traces_());
%             uimenu(mh_view,'Label', 'View all [R]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'r'));
%             uimenu(mh_view,'Label', '[Z]oom selected', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'z'));
%             uimenu(mh_view,'Label', '[W]aveform (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'w'));
%             uimenu(mh_view,'Label', '[N]umbers (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'n'));
%             uimenu(mh_view,'Label', 'Show raw waveform', 'Callback', @(h,e)raw_waveform_(h), ...
%             'Checked', ifeq_(get_(obj.hCfg, 'showRaw'), 'on', 'off'));
%             %uimenu(mh_view,'Label', 'Threshold by sites', 'Callback', @(h,e)keyPressFcn_thresh_(hFig, 'n'));
%             % uimenu(mh_view,'Label', '.prm file', 'Callback', @edit_prm_);
%             uimenu(mh_view,'Label', 'Reset window positions', 'Callback', @reset_position_);
% 
%             mh_proj = uimenu(hFig,'Label','Projection');
%             uimenu(mh_proj, 'Label', 'vpp', 'Callback', @(h,e) proj_view_(h), ...
%             'Checked', if_on_off_(obj.hCfg.vcFet_show, {'vpp', 'vmin'}));
%             uimenu(mh_proj, 'Label', 'pca', 'Callback', @(h,e) proj_view_(h), ...
%             'Checked', if_on_off_(obj.hCfg.vcFet_show, {'pca'}));
%             uimenu(mh_proj, 'Label', 'ppca', 'Callback', @(h,e) proj_view_(h), ...
%             'Checked', if_on_off_(obj.hCfg.vcFet_show, {'ppca', 'private pca'}));
%             % uimenu(mh_proj, 'Label', 'cov', 'Callback', @(h,e)proj_view_(h), ...
%             %     'Checked', if_on_off_(P.vcFet_show, {'cov', 'spacetime'}));
% 
%             mh_plot = uimenu(hFig,'Label','Plot');
%             uimenu(mh_plot, 'Label', 'All unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_);
%             uimenu(mh_plot, 'Label', 'Selected unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_(1));
% 
            obj.hMenus('hMenuInfo') = uimenu(hFig, 'Label', '', 'Tag', 'hMenuInfo');
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'Annotate unit', 'Callback', @unit_annotate_);
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'single', 'Callback', @(h,e)unit_annotate_(h,e,'single'));
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'multi', 'Callback', @(h,e)unit_annotate_(h,e,'multi'));
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'noise', 'Callback', @(h,e)unit_annotate_(h,e,'noise'));
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'clear annotation', 'Callback', @(h,e)unit_annotate_(h,e,''));
%             uimenu(obj.hMenus('hMenuInfo'), 'Label', 'equal to', 'Callback', @(h,e)unit_annotate_(h,e,'=%d'));
% 
%             mh_help = uimenu(hFig,'Label','Help');
%             uimenu(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);

            drawnow;
            hFig.outerPosition = outerPosition;
        end

        function closeFigures(obj)
            %CLOSEFIGURES Close all open figures
            if obj.hasFig('hFigWav') && obj.hFigs('hFigWav').isReady
                hFigWav = obj.hFigs('hFigWav');
                remove(obj.hFigs, 'hFigWav'); % to prevent infinite recursion!
                hFigWav.close(); % calls killFigWav
            end

            try
                obj.figApply(@(hFig) hFig.close());
            catch ME
            end

            obj.hFigs = containers.Map();
        end

        function deleteClusters(obj, clusters)
            if numel(obj.selected) > 1 && nargin < 2
                return;
            elseif nargin == 2
                deleteMe = clusters;
            else
                deleteMe = obj.selected(1);
            end

            success = obj.hClust.deleteClusters(deleteMe);
            if success
                % save the new clustering
                deleted = strjoin(arrayfun(@num2str, deleteMe, 'UniformOutput', false), ', ');
                commitMsg = sprintf('%s: delete %s', datestr(now, 31), deleted);
                obj.hClust.commit(commitMsg);

                % replot
                doPlotFigWav(obj.hFigs('hFigWav'), obj.hClust, obj.hCfg, obj.maxAmp);
                doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);
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

        function killFigWav(obj, hObject, ~)
            %KILLFIGWAV Destroy the main figure, close all other figures
            if ~obj.isEnding
                obj.endSession(); % we'll be back
            end

            delete(hObject);
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
                commitMsg = sprintf('%s: merge %d into %d', datestr(now, 31), jCluster, iCluster);
                obj.hClust.commit(commitMsg);

                % replot
                doPlotFigWav(obj.hFigs('hFigWav'), obj.hClust, obj.hCfg, obj.maxAmp);
                doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);
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

            % plot rho-delta figure
            if obj.hasFig('hFigRD')
                doPlotFigRD(obj.hFigs('hFigRD'), obj.hClust, obj.hCfg);
            end

            % plot sim score figure
            if obj.hasFig('hFigSim')
                hFigSim = doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);

                % set key and mouse handles
                hFigSim.hFunKey = @obj.keyPressFigSim;
                hFigSim.setMouseable(@obj.mouseClickFigSim);
            end

%                 case 'u', update_FigCor_(S0);
%                 case 'p' %PSTH plot
%                 if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
%                 plot_raster_(S0, 1);

            % plot feature projection
            if obj.hasFig('hFigProj')
                % set key and mouse handles
                hFigProj = obj.hFigs('hFigProj');
                hFigProj.hFunKey = @obj.keyPressFigProj;
                hFigProj.setMouseable(); % nothing special here

                % obj.hFigs('hFigProj') = hFigProj;
            end

            % plot amplitude vs. time
            if obj.hasFig('hFigTime')
                hFigTime = doPlotFigTime(obj.hFigs('hFigTime'), obj.hClust, obj.hCfg, obj.selected, obj.maxAmp);

                % set key and mouse handles
                hFigTime.hFunKey = @obj.keyPressFigTime;
                hFigTime.setMouseable(); % nothing special here

                % obj.hFigs('hFigTime') = hFigTime;
            end

            % plot main waveform view
            if obj.hasFig('hFigWav')
                hFigWav = doPlotFigWav(obj.hFigs('hFigWav'), obj.hClust, obj.hCfg, obj.maxAmp);

                % set key and mouse handles
                hFigWav.hFunKey = @obj.keyPressFigWav;
                hFigWav.setMouseable(@obj.mouseClickFigWav);

                % make this guy the key log
                hFigWav.figSet('CloseRequestFcn', @obj.killFigWav);
                obj.addMenu(hFigWav);

                % obj.hFigs('hFigWav') = hFigWav;
            end

            % select first cluster
            obj.updateSelect(1);
        end

        function updateSelect(obj, iClusters)
            %UPDATESELECT Select a cluster or pair of clusters across all views
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
            obj.currentSite = obj.hClust.clusterSites(obj.selected(1));

            % update plots
            obj.updateCursorFigWav();
            obj.updateCursorFigSim();
            obj.updateFigCorr();
            obj.updateFigHist();
            obj.updateFigISI();
            obj.updateFigMap();
            obj.updateFigPos();
            obj.updateFigProj();
            obj.updateFigTime();

            % autoscale figProj and figTime
            autoScaleProjTime(obj.hClust, obj.hFigs('hFigProj'), obj.hFigs('hFigTime'), obj.selected);

            % update menu entry to indicate selected clusters
            if numel(obj.selected) > 1 && obj.hasMenu('hMenuInfo')
                menuLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', obj.selected(1), ...
                                    obj.hClust.clusterNotes{obj.selected(1)}, obj.selected(2), ...
                                    obj.hClust.clusterNotes{obj.selected(2)});
                set(obj.hMenus('hMenuInfo'), 'Label', menuLabel);
            elseif obj.hasMenu('hMenuInfo')
                menuLabel = sprintf('Unit %d "%s"', obj.selected(1), obj.hClust.clusterNotes{obj.selected(1)});
                set(obj.hMenus('hMenuInfo'), 'Label', menuLabel);
            end

            % zoom in on waveform view
            %obj.keyPressFigWav([], struct('Key', 'z'));
        end

        function spawnFigures(obj)
            %SPAWNFIGURES Create new figures
            obj.hFigs = doSpawnFigures(obj.hCfg);
        end
    end

    %% USER METHODS
    methods
        function allToForeground(obj)
            %ALLTOFOREGROUND Bring all figures to the foreground
            obj.figApply(@(hFig) hFig.toForeground());
        end

        function beginSession(obj, hClust)
            %BEGINSESSION Start curating clusters
            if ~isempty(obj.hFigs) % session already running
                return;
            end
            obj.isEnding = false;
            obj.cRes = struct('hClust', hClust);
            obj.selected = 1;
            obj.currentSite = hClust.clusterSites(1);
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
        % hClust
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = [];
            end
        end

%         % selected
%         function set.selected(obj, se)
%             obj.figApply(@(hFig) setfield(hFig.figData.selected, se));
%             obj.selected = se;
%         end
    end
end

