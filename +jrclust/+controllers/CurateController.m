classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters

    properties (Access=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hCfg;           % Config object
        hClust;         % Clustering object
    end

    properties (AbortSet, Access=private, Hidden, Transient, SetObservable)
        hFigs;          % containers.Map of Figure objects
        selected;       % selected clusters, in order of selection
    end

    %% LIFECYCLE
    methods
        function obj = CurateController(hCfg)
            obj.hCfg = hCfg;
            obj.hFigs = containers.Map();
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
        function keyPressFigSim(obj, hFigSim, hEvent)
            %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
            disp(hEvent);
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

        function keyPressFigWav(obj, hObject, hEvent)
            %KEYPRESSFIGWAV Handle callbacks for keys pressed in main view
            hFigWav = obj.hFigs('hFigWav');
            figData = hFigWav.figGet('UserData');

            nSites = numel(obj.hCfg.siteMap);

            switch hEvent.Key
                case {'uparrow', 'downarrow'} % rescale
                    disp('rescale');
                    %rescale_FigWav_(hEvent, S0, hCfg);
                    %clu_info_(S0); %update figpos

                case 'leftarrow' % select previous cluster
                    disp('previous cluster');

                case 'rightarrow' % select next cluster
                    disp('next cluster');

                case 'home' % select first cluster
                    disp('first cluster');

                case 'end' % select last cluster
                    disp('last cluster');
%                     if strcmpi(hEvent.Key, 'home')
%                         S0.iCluCopy = 1;
%                     elseif strcmpi(hEvent.Key, 'end')
%                         S0.iCluCopy = hClust.nClusters;
%                     elseif ~key_modifier_(hEvent, 'shift');
%                         if strcmpi(hEvent.Key, 'leftarrow')
%                             if S0.iCluCopy == 1, return; end
%                             S0.iCluCopy = S0.iCluCopy - 1;
%                         else
%                             if S0.iCluCopy == hClust.nClusters, return; end
%                             S0.iCluCopy = S0.iCluCopy + 1;
%                         end
%                     else
%                         if isempty(S0.iCluPaste)
%                             S0.iCluPaste = S0.iCluCopy;
%                         end
%                         if strcmpi(hEvent.Key, 'leftarrow')
%                             if S0.iCluPaste == 1, return; end
%                             S0.iCluPaste = S0.iCluPaste - 1;
%                         else
%                             if S0.iCluPaste == hClust.nClusters, return; end
%                             S0.iCluPaste = S0.iCluPaste + 1;
%                         end
%                     end
%                     S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0); %select first clu
%                     if strcmpi(hEvent.Key, 'home') || strcmpi(hEvent.Key, 'end') %'z' to recenter
%                         S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0);
%                     end
                case 'm' % merge clusters
                    % S0 = ui_merge_(S0); % merge clusters
                    disp('merge');

                case 'space' % select most similar to currently selected
                    disp('select most similar');
%                     mrWavCor = hClust.mrWavCor;
%                     mrWavCor(S0.iCluCopy,S0.iCluCopy) = -inf;
%                     [~,S0.iCluPaste] = max(mrWavCor(:,S0.iCluCopy));
%                     set(0, 'UserData', S0);
%                     button_CluWav_simulate_([], S0.iCluPaste);

                case 's' % split
                    % auto_split_(1, S0);
                    disp('split');

                case 'r' %reset view
                    %figure_wait_(1);
                    %axis_([0, S0.S_clu.nClusters + 1, 0, numel(hCfg.viSite2Chan) + 1]);
                    %figure_wait_(0);
                    disp('reset');

                case {'d', 'backspace', 'delete'}
                    % S0 = ui_delete_(S0);
                    disp('delete');

                case 'z' % zoom
                    if isempty(obj.selected)
                        obj.selected = 1;
                    end
                    iSite = obj.hClust.clusterSites(obj.selected(1));
                    hFigWav.setWindow(obj.selected(1) + [-1, 1]*6, iSite + [-1, 1]*(obj.hCfg.maxSite*2+1), [0 obj.hClust.nClusters+1], [0 nSites+1]);

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
%                 case 't', plot_FigTime_(S0); % time view
%                 case 'j', plot_FigProj_(S0); %projection view
%                 case 'n'
%                 fText = get_set_(figData, 'fText', get_set_(hCfg, 'fText', 1));
%                 setFigWavXTicks(figData, hClust, ~fText);
%                 case 'i', plot_FigHist_(S0); %ISI histogram
%                 case 'e', plot_FigMap_(S0);
%                 case 'u', update_FigCor_(S0);
%                 case 'p' %PSTH plot
%                 if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
%                 plot_raster_(S0, 1);
                otherwise
                    hFigWav.wait(false); %stop waiting
            end

            hFigWav.toForeground(); %change the focus back to the current object
        end

        function mouseClickFigSim(obj, xyPos, clickType)
            %MOUSECLICKFIGSIM Handle callbacks for mouse clicks in sim view
            xyPos = floor(xyPos);
            if strcmp(clickType, 'normal') % left click
                obj.selected = xyPos(1); % first selected cluster is x position

                if diff(xyPos) ~= 0
                    obj.selected(2) = xyPos(2);
                end

                disp(obj.selected);

                % S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
                % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
            end
        end

        function mouseClickFigWav(obj, xyPos, clickType)
            %MOUSECLICKFIGWAV Handle callbacks for mouse clicks in sim view
            iCluster = floor(xyPos(1)); % floor of x position
            if iCluster < 1 || iCluster > obj.hClust.nClusters
                return;
            end

            if strcmp(clickType, 'normal')  % left click, select primary cluster
                obj.selected = iCluster;
                obj.updateCursorFigWav();
            elseif strcmp(clickType, 'alt') && iCluster ~= obj.selected(1) % right click, select secondary cluster
                obj.selected(2) = iCluster;
                obj.updateCursorFigWav();
            else                            % middle click, ignore
                disp(clickType);
                return;
            end

            hFig = obj.hFigs('hFigWav');
            hFig.wait(true);
            % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z'
            % auto_scale_proj_time_(S0);
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
                colorMap = [1 0 0]; % red
            elseif strcmp(plotKey, 'selected1')
                colorMap = [0 0 0]; % black
            else
                return; % maybe some day we support selected3, 4, ...
            end

            if ~hFigWav.hasPlot(plotKey)
                hFigWav.addPlot(plotKey, nan, nan, 'Color', colorMap, 'LineWidth', 2);
            end

            if obj.hCfg.fWav_raw_show
                meanWf = obj.hClust.meanWfGlobalRaw(:, :, iCluster);
            else
                meanWf = obj.hClust.meanWfGlobal(:, :, iCluster);
            end

            % this is a hack until we hunt down and destroy all calls to
            % `multiplot` which are actually rescales (nargin <= 2)
            userData = hFigWav.plotGet(plotKey, 'UserData');
            if isempty(userData)
                userData = struct('maxAmp', obj.hCfg.maxAmp);
            elseif ~isfield(userData, 'maxAmp')
                userData.maxAmp = obj.hCfg.maxAmp;
            end 

            hFigWav.multiplot(plotKey, userData.maxAmp, getXRange(iCluster, obj.hCfg), meanWf);
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
                colorH = [0 0 0]; % black
            else
                colorH = [1 0 0]; % red
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
                obj.plotSelectedWaveforms(obj.selected(2), 'selected2');
            end

            obj.updateCursorFigSim();
            obj.updateFigCorr();
        end

        function updateFigCorr(obj)
            %UPDATEFIGCORR Plot cross correlation
            if isempty(obj.selected) || ~obj.hasFig('hFigCorr')
                return;
            end

            hFigCorr = doPlotFigCorr(obj.hFigs('hFigCorr'), obj.hClust, obj.hCfg, obj.selected);
            obj.hFigs('hFigCorr') = hFigCorr;
        end
        %                 case 'i', plot_FigHist_(S0); %ISI histogram

        function updateFigHist(obj)
            %UPDATEFIGHIST Plot ISI histogram
            if isempty(obj.selected) || ~obj.hasFig('hFigHist')
                return;
            end
            
            hFigHist = doPlotFigHist(obj.hFigs('hFigHist'), obj.hClust, obj.hCfg, obj.selected);
            obj.hFigs('hFigHist') = hFigHist;
        end

        function updateFigISI(obj)
            %UPDATEFIGISI Plot return map
            if isempty(obj.selected) || ~obj.hasFig('hFigISI')
                return;
            end

            hFigISI = doPlotFigISI(obj.hFigs('hFigISI'), obj.hClust, obj.hCfg, obj.selected);
            obj.hFigs('hFigISI') = hFigISI;
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        function addMenu(obj, hFig)
            drawnow;

            outerPosition = hFig.outerPosition;
            hFig.figSet('MenuBar','None');

            hMFile = hFig.uimenu('Label', 'File');
            uimenu(hMFile, 'Label', 'Save', 'Callback', @obj.saveFiles); % save_manual_
            uimenu(hMFile, 'Label', 'Save figures as .fig', 'Callback', @(h, e) obj.saveFigures('.fig')); % save_figures_
            uimenu(hMFile, 'Label', 'Save figures as .png', 'Callback', @(h, e) obj.saveFigures('.png')); % save_figures_
%             uimenu(mhFile, 'Label', 'Describe', 'Callback', @(h,e) msgbox_(describe_()), 'Separator', 'on');
%             uimenu(mhFile, 'Label', 'Edit prm file', 'Callback', @edit_prm_);
%             uimenu(mhFile, 'Label', 'Reload prm file', 'Callback', @reload_prm_);
%             uimenu(mhFile, 'Label', 'Export units to csv', 'Callback', @export_csv_, 'Separator', 'on');
%             uimenu(mhFile, 'Label', 'Export unit qualities to csv', 'Callback', @(h,e)export_quality_);
%             uimenu(mhFile, 'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
%             uimenu(mhFile, 'Label', 'Export selected mean unit waveforms', 'Callback', @(h,e)export_mrWav_clu_);
%             uimenu(mhFile, 'Label', 'Export all waveforms from the selected unit', 'Callback', @(h,e)export_tnWav_spk_);
%             uimenu(mhFile, 'Label', 'Export firing rate for all units', 'Callback', @(h,e)export_rate_);
            uimenu(hMFile, 'Label', 'Exit', 'Callback', @(h, e) obj.endSession, 'Separator', 'on', 'Accelerator', 'Q');

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
%             'Checked', ifeq_(get_(obj.hCfg, 'fWav_raw_show'), 'on', 'off'));
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
%             mh_info = uimenu(hFig,'Label','','Tag', 'mh_info');
%             uimenu(mh_info, 'Label', 'Annotate unit', 'Callback', @unit_annotate_);
%             uimenu(mh_info, 'Label', 'single', 'Callback', @(h,e)unit_annotate_(h,e,'single'));
%             uimenu(mh_info, 'Label', 'multi', 'Callback', @(h,e)unit_annotate_(h,e,'multi'));
%             uimenu(mh_info, 'Label', 'noise', 'Callback', @(h,e)unit_annotate_(h,e,'noise'));
%             uimenu(mh_info, 'Label', 'clear annotation', 'Callback', @(h,e)unit_annotate_(h,e,''));
%             uimenu(mh_info, 'Label', 'equal to', 'Callback', @(h,e)unit_annotate_(h,e,'=%d'));
% 
%             mh_help = uimenu(hFig,'Label','Help');
%             uimenu(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);

            drawnow;
            hFig.outerPosition = outerPosition;
        end

        function closeFigures(obj)
            %CLOSEFIGURES Close all open figures
            if obj.hasFig('hFigWav')
                hFigWav = obj.hFigs('hFigWav');
                hFigWav.close(); % calls killFigWav
            else
                try
                    obj.figApply(@(hFig) hFig.close());
                catch ME
                end

                obj.hFigs = containers.Map();
            end
        end

        function res = figApply(obj, hFun)
            %FIGAPPLY Apply a function to all figures
            if nargout == 0 % "Too many output arguments"
                cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs));
            else
                res = cellfun(@(k) hFun(obj.hFigs(k)), keys(obj.hFigs));
            end
        end

        function killFigWav(obj, hObject, ~)
            %KILLFIGWAV Destroy the main figure, close all other figures
            if obj.hasFig('hFigWav')
                remove(obj.hFigs, 'hFigWav'); % to prevent infinite recursion!
            end

            delete(hObject);
            obj.closeFigures();
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
                obj.hFigs('hFigRD') = doPlotFigRD(obj.hFigs('hFigRD'), obj.hClust, obj.hCfg);
            end

            % plot sim score figure
            if obj.hasFig('hFigSim')
                hFigSim = doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);

                % set key and mouse handles
                hFigSim.hFunKey = @obj.keyPressFigSim;
                hFigSim.setMouseable(@obj.mouseClickFigSim);

                obj.hFigs('hFigSim') = hFigSim;
            end

            % plot main waveform view
            if obj.hasFig('hFigWav')
                hFigWav = doPlotFigWav(obj.hFigs('hFigWav'), obj.hClust, obj.hCfg);

                % set key and mouse handles
                hFigWav.hFunKey = @obj.keyPressFigWav;
                hFigWav.setMouseable(@obj.mouseClickFigWav);

                % make this guy the key log
                hFigWav.figSet('CloseRequestFcn', @obj.killFigWav);
                obj.addMenu(hFigWav);

                obj.hFigs('hFigWav') = hFigWav;
            end

            % plot feature projection view
%             if obj.hasFig('hFigProj')
%                 
%             end

%                 case 't', plot_FigTime_(S0); % time view
%                 case 'j', plot_FigProj_(S0); %projection view
%                 case 'e', plot_FigMap_(S0);
%                 case 'u', update_FigCor_(S0);
%                 case 'p' %PSTH plot
%                 if isempty(hCfg.vcFile_trial), msgbox_('''vcFile_trial'' not set. Reload .prm file after setting (under "File menu")'); return; end
%                 plot_raster_(S0, 1);

            % plot correlation view
            obj.updateFigCorr();

            % plot ISI histogram
            obj.updateFigHist();

            % plot return map (ISI view)
            obj.updateFigISI();

            % zoom in on waveform view
            obj.keyPressFigWav([], struct('Key', 'z'));
        end

        function spawnFigures(obj)
            %SPAWNFIGURES Create new figures
            obj.hFigs = doSpawnFigures(obj.hCfg);
        end
    end

    %% USER METHODS
    methods
        function beginSession(obj, hClust)
            %BEGINSESSION Start curating clusters
            obj.cRes = struct('hClust', hClust);
            obj.selected = 1;
            obj.plotAllFigures();
        end

        function res = endSession(obj)
            %ENDSESSION Finish curating and return results
            if nargout == 0
                obj.saveFiles(); % ask
            end
            obj.closeFigures();

            obj.selected = [];
            res = obj.cRes;
            obj.cRes = [];
        end

        function hf = hasFig(obj, figKey)
            %HASFIG Return true if we have a figure by key
            hf = ischar(figKey) && isKey(obj.hFigs, figKey);
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

