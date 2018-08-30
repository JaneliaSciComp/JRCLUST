%--------------------------------------------------------------------------
function keyPressFigProj(hFig, event)
    S0 = get(0, 'UserData');
    P = S0.P;

    [hFig, figData] = getCachedFig('FigProj');

    plotDataFG = get(figData.hPlotFG, 'UserData');
    sitesOfInterest = plotDataFG.sitesOfInterest;

    figure_wait_(1);

    switch lower(event.Key)
        case {'uparrow', 'downarrow'}
            rescaleFigProj(event, hFig, figData, S0);

        case {'leftarrow', 'rightarrow'} % change channels
            fPlot = 0;
            if strcmpi(event.Key, 'leftarrow')
                if min(figData.sitesOfInterest) > 1
                    figData.sitesOfInterest = figData.sitesOfInterest - 1;
                    fPlot = 1;
                end
            else
                if max(figData.sitesOfInterest) < max(P.chanMap)
                    figData.sitesOfInterest = figData.sitesOfInterest + 1;
                    fPlot = 1;
                end
            end

            if fPlot
                % need to refresh features for this set of sites
                set(hFig, 'UserData', figData);
                S0.P.sitesOfInterest = figData.sitesOfInterest;
                plotFigProj(S0);
            end

        case 'r' % reset view
            nSites = min(getOr(P, 'nSitesFigProj', 5), numel(sitesOfInterest));
            axis_([0 nSites 0 nSites]);
            plotFigProj(S0);

        case 's' % split
            figure_wait_(0);
            if ~isempty(S0.secondarySelectedCluster)
                msgbox_('Select one cluster to split');
                return;
            end

            clearPlots(figData.hPlotFG2); % replacement for update_plot2_proj_ in select_polygon_
            plotDataFG = select_polygon_(figData.hPlotFG);

            if ~isempty(plotDataFG)
                [fSplit, vlIn] = plot_split_(plotDataFG);
                if fSplit
                    S_clu = splitCluster(S0.primarySelectedCluster, vlIn);
                else
                    clearPlots(figData.hPlotFG2); % replacement for update_plot2_proj_
                    % update_plot2_proj_();
                end
            end

        case 'm'
            manualMerge(S0);

        case 'f'
            if getOr(P, 'fImportKilosort')
                if strcmpi(P.displayFeature, 'vpp')
                    P.displayFeature = 'kilosort';
                else
                    P.displayFeature = 'vpp';
                    figData.maxAmp = P.maxAmp;
                    set(hFig, 'UserData', figData);
                end
                S0.P = P;
                set(0, 'UserData', S0);

                plotFigProj(S0);
            else
                disp('keyPressFigProj: ''f'': not implemented yet');
            end

        case 'p' % toggle PCivPCj
            if getOr(P, 'fImportKilosort') && strcmpi(P.displayFeature, 'kilosort')
                if all(S0.pcPair == [1 2])
                    S0.pcPair = [1 3];
                elseif all(S0.pcPair == [1 3])
                    S0.pcPair = [2 3];
                else
                    S0.pcPair = [1 2];
                end
                set(0, 'UserData', S0);

                plotFigProj(S0);
            end

        case 'b' %background spikes
            toggleVisible_(figData.hPlotBG);

        case 'h' %help
            msgbox_(figData.csHelp, 1);
    end % switch

    figure_wait_(0);
end % function
