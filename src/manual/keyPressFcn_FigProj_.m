%--------------------------------------------------------------------------
function keyPressFcn_FigProj_(hFig, event)
    S0 = get(0, 'UserData');
    [P, S_clu] = deal(S0.P, S0.S_clu);
    [hFig, S_fig] = get_fig_cache_('FigProj');
    % nSites = numel(P.viSite2Chan);
    S_plot1 = get(S_fig.hPlot1, 'UserData');
    viSites_show = S_plot1.viSites_show;
    % nSites = numel(viSites_show);
    % set(hObject, 'Pointer', 'watch');
    figure_wait_(1);
    switch lower(event.Key)
        case {'uparrow', 'downarrow'}
            rescaleFigProj(event, hFig, S_fig, S0);

        case {'leftarrow', 'rightarrow'} % change channels
            fPlot = 0;
            if strcmpi(event.Key, 'leftarrow')
                if min(S_fig.viSites_show)>1
                    S_fig.viSites_show=S_fig.viSites_show-1;
                    fPlot = 1;
                end
            else
                if max(S_fig.viSites_show) < max(P.viSite2Chan)
                    S_fig.viSites_show=S_fig.viSites_show+1;
                    fPlot = 1;
                end
            end
            if fPlot
                set(hFig, 'UserData', S_fig);
                S0.P.viSites_show = S_fig.viSites_show;
                plot_FigProj_(S0);
            end

        case 'r' %reset view
            axis_([0 numel(viSites_show) 0 numel(viSites_show)]);

        case 's' %split
            figure_wait_(0);
            if ~isempty(S0.iCluPaste)
                msgbox_('Select one cluster to split'); return;
            end
            S_plot1 = select_polygon_(S_fig.hPlot1);
            if ~isempty(S_plot1)
                [fSplit, vlIn] = plot_split_(S_plot1);
                if fSplit
                    S_clu = split_clu_(S0.iCluCopy, vlIn);
                else
                    update_plot2_proj_();
                    %                 delete_multi_(S_plot1.hPoly);
                end
            end

        case 'm'
            ui_merge_(S0);

        case 'f'
            if get_set_(P, 'fImportKsort', 0) && strcmpi(P.vcFet_show, 'vpp')
                P.vcFet_show = 'kilosort';
            elseif strcmpi(P.vcFet_show, 'vpp')
                P.vcFet_show = P.vcFet;
            else
                P.vcFet_show = 'vpp';
            end

            set0_(P);
            try
                plot_FigProj_();
            catch
                P.vcFet_show = 'vpp';
                set0_(P);
                plot_FigProj_();
            end
            
        case 'p' % toggle PCi v. PCj
            if get_set_(P, 'fImportKsort', 0) && strcmpi(P.vcFet_show, 'kilosort')
                pcPair = get_set_(S0, 'pcPair', [1 2]);

                if all(pcPair == [1 2])
                    S0.pcPair = [1 3];
                elseif all(pcPair == [1 3])
                    S0.pcPair = [2 3];
                else
                    S0.pcPair = [1 2];
                end

                set(0, 'UserData', S0);
                
                plot_FigProj_(S0);
            end

        case 'b' % background spikes
            toggleVisible_(S_fig.hPlot0);

        case 'h' % help
            msgbox_(S_fig.csHelp, 1);
    end % switch

    % drawnow;
    figure_wait_(0);
end %func
