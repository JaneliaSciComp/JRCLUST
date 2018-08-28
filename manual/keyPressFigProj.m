%--------------------------------------------------------------------------
function keyPressFigProj(hFig, event)
    S0 = get(0, 'UserData');
    P = S0.P;
    S_clu = S0.S_clu;

    [hFig, S_fig] = getCachedFig('FigProj');

    S_plot1 = get(S_fig.hPlotFG, 'UserData');
    sitesOfInterest = S_plot1.sitesOfInterest;

    figure_wait_(1);
    switch lower(event.Key)
        case {'uparrow', 'downarrow'}
            rescale_FigProj_(event, hFig, S_fig, S0);

        case {'leftarrow', 'rightarrow'} % change channels
            fPlot = 0;
            if strcmpi(event.Key, 'leftarrow')
                if min(S_fig.sitesOfInterest) > 1
                    S_fig.sitesOfInterest = S_fig.sitesOfInterest - 1;
                    fPlot = 1;
                end
            else
                if max(S_fig.sitesOfInterest) < max(P.chanMap)
                    S_fig.sitesOfInterest = S_fig.sitesOfInterest + 1;
                    fPlot = 1;
                end
            end

            if fPlot
                set(hFig, 'UserData', S_fig);
                S0.P.sitesOfInterest = S_fig.sitesOfInterest;
                plotFigProj(S0);
            end

        case 'r' % reset view
            nSites = min(getOr(P, 'nSitesFigProj', 5), numel(sitesOfInterest));
            axis_([0 nSites 0 nSites]);
            plotFigProj(S0);

        case 's' % split
            figure_wait_(0);
            if ~isempty(S0.secondarySelectedCluster)
                msgbox_('Select one cluster to split'); return;
            end

            S_plot1 = select_polygon_(S_fig.hPlotFG);

            if ~isempty(S_plot1)
                [fSplit, vlIn] = plot_split_(S_plot1);
                if fSplit
                    S_clu = splitCluster(S0.primarySelectedCluster, vlIn);
                else
                    update_plot2_proj_();
                end
            end

        case 'm'
            ui_merge_(S0);

        case 'f'
            disp('keyPressFigProj: ''f'': not implemented yet');
            % if strcmpi(P.displayFeature, 'vpp')
            %     P.displayFeature = P.feature;
            % else
            %     P.displayFeature = 'vpp';
            % end
            %
            % S0.P = P;
            % set(0, 'UserData', S0);
            %
            % plotFigProj();

        case 'p' % toggle PCivPCj
            if getOr(P, 'fImportKilosort')
                if S0.kspc == [1 2]
                    S0.kspc = [1 3];
                elseif S0.kspc == [1 3]
                    S0.kspc = [2 3];
                else
                    S0.kspc = [1 2];
                end
                set(0, 'UserData', S0);

                plotFigProj(S0);
            end

        case 'b' %background spikes
            toggleVisible_(S_fig.hPlotBG);

        case 'h' %help
            msgbox_(S_fig.csHelp, 1);
    end %switch
    % drawnow;
    figure_wait_(0);
end %func
