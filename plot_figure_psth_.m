%--------------------------------------------------------------------------
function plot_figure_psth_(hFig, iClu, crTime_trial, S_clu, P)
    S_fig = get(hFig, 'UserData');
    [vhAx1, vhAx2, vcColor] = deal(S_fig.vhAx1, S_fig.vhAx2, S_fig.vcColor);
    for iStim = 1:numel(vhAx1)
        cla(vhAx1(iStim));
        cla(vhAx2(iStim));
        vrTime_trial = crTime_trial{iStim}; %(:,1);
        nTrials = numel(vrTime_trial);
        viTime_clu1 = S_clu_time_(S_clu, iClu);
        plot_raster_clu_(viTime_clu1, vrTime_trial, P, vhAx1(iStim));
        plot_psth_clu_(viTime_clu1, vrTime_trial, P, vhAx2(iStim), vcColor);
        title(vhAx2(iStim), sprintf('Cluster %d; %d trials', iClu, nTrials));
    end
    %     offset = offset + nTrials;
    if numel(vhAx1)>2
        set(vhAx1(2:end),'xticklabel',{});
        for ax = vhAx1(2:end)
            xlabel(ax, '')
        end
    end % end
end % function
