classdef CurateController < handle
    %CURATECONTROLLER Interface for manually curating sorted clusters

    properties (Access=private, Hidden, SetObservable)
        cRes;           % curate results struct, returned at endSession
        hCfg;           % Config object
        hClust;         % Clustering object
    end

    properties (Access=private, Hidden, Transient)
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
        function keyPressFigSim(obj, ~, hEvent)
            %KEYPRESSFIGSIM Handle callbacks for keys pressed in sim view
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
            %disp(hObject);
            %disp(hEvent);
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
            disp(clickType);
            if strcmp(clickType, 'normal')  % select primary cluster
                obj.updateCursor(floor(xyPos(1)), false)
            elseif strcmp(clickType, 'alt') % select secondary cluster
                obj.updateCursor(floor(xyPos(1)), true)
            else
                return;
            end

            hFig = obj.hFigs('hFigWav');
            hFig.figSet('Pointer', 'watch');
            % S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j','t','c','i','v','e','f'}, S0); %'z'
            % auto_scale_proj_time_(S0);
            % plot_raster_(S0);
            hFig.figSet('Pointer', 'arrow');
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
            hKeys = keys(obj.hFigs);
            if any(strcmp(hKeys, 'hFigWav'))
                hFigWav = obj.hFigs('hFigWav');
                hFigWav.close(); % calls killFigWav, removes itself, and calls this function
            else
                for iKey = 1:numel(hKeys)
                    hFig = obj.hFigs(hKeys{iKey});
                    try
                        hFig.close();
                    catch ME
                    end
                end

                obj.hFigs = containers.Map();
            end
        end

        function killFigWav(obj, hObject, hEvent)
            hKeys = keys(obj.hFigs);
            if any(strcmp(hKeys, 'hFigWav'))
                remove(obj.hFigs, 'hFigWav');
            end

            delete(hObject);
            obj.closeFigures();
        end

        function plotAllFigures(obj)
            %PLOTALLFIGURES Plot all figures
            if isempty(obj.hFigs)
                obj.spawnFigures();
            end
            hKeys = keys(obj.hFigs);

            % plot rho-delta figure
            if any(strcmp(hKeys, 'hFigRD'))
                obj.hFigs('hFigRD') = doPlotFigRD(obj.hFigs('hFigRD'), obj.hClust, obj.hCfg);
            end

            % plot sim score figure
            if any(strcmp(hKeys, 'hFigSim'))
                hFigSim = doPlotFigSim(obj.hFigs('hFigSim'), obj.hClust, obj.hCfg);

                % set key and mouse handles
                hFigSim.hFunKey = @obj.keyPressFigSim;
                hFigSim.setMouseable(@obj.mouseClickFigSim);

                obj.hFigs('hFigSim') = hFigSim;
            end

            % plot main waveform view
            if any(strcmp(hKeys, 'hFigWav'))
                hFigWav = doPlotFigWav(obj.hFigs('hFigWav'), obj.hClust, obj.hCfg);

                % set key and mouse handles
                hFigWav.hFunKey = @obj.keyPressFigWav;
                hFigWav.setMouseable(@obj.mouseClickFigWav);

                % make this guy the main view
                hFigWav.figSet('CloseRequestFcn', @obj.killFigWav);
                obj.addMenu(hFigWav);

                obj.hFigs('hFigWav') = hFigWav;
            end
        end

        function spawnFigures(obj)
            %SPAWNFIGURES Create new figures
            obj.hFigs = doSpawnFigures(obj.hCfg);
        end

        function updateCursor(obj, iCluster, fPaste)
            if isempty(iCluster)
                return;
            end
%             if ~isfield(S0, 'hCopy'), S0.hCopy = []; end
%             if ~isfield(S0, 'hPaste'), S0.hPaste = []; end

            if ~fPaste
                selected_ = iCluster;
                if selected_ < 1 || selected_ > obj.hClust.nClusters
                    return;
                end
                % update_plot_(S0.hPaste, nan, nan); %hide paste
                obj.selected = selected_;
                % [S0.iCluCopy, S0.hCopy] = plot_tmrWav_clu_(S0, iCluCopy, S0.hCopy, [0 0 0]);
            else
                selected_ = obj.selected;
                selected_(2) = iCluster;
                if selected_(2) < 1 || selected_(2) > obj.hClust.nClusters || numel(unique(selected_)) ~= numel(selected_)
                    return;
                end
                obj.selected = selected_;
                % [S0.iCluPaste, S0.hPaste] = plot_tmrWav_clu_(S0, iCluPaste, S0.hPaste, [1 0 0]);
            end
            % set(hFig, 'UserData', S_fig);
            % cursor_FigWavCor_(S0);
            % if nargout==0, set(0, 'UserData', S0); end
        end
    end

    %% USER METHODS
    methods
        function beginSession(obj, hClust)
            %BEGINSESSION Start curating clusters
            obj.cRes = struct('hClust', hClust);
            obj.plotAllFigures();
        end

        function res = endSession(obj)
            %ENDSESSION Finish curating and return results
            if nargout == 0
                obj.saveFiles(); % ask
            end
            obj.closeFigures();

            res = obj.cRes;
            obj.cRes = [];
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
        function hc = get.hClust(obj)
            if ~isempty(obj.cRes) && isfield(obj.cRes, 'hClust')
                hc = obj.cRes.hClust;
            else
                hc = [];
            end
        end
    end
end

