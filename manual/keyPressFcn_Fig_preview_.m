%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function S0 = keyPressFcn_Fig_preview_(hFig, event, S0)
    if nargin<3, S0 = get(0, 'UserData'); end
    P = get_(S0, 'P');
    S_fig = get(hFig, 'UserData');
    factor = 1 + 3 * keyModifier(event, 'shift');
    nSites = numel(P.chanMap);

    switch lower(event.Key)
        case {'uparrow', 'downarrow'}
        S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, ...
        S_fig.hPlot_traces, S_fig.hPlot_traces_bad, S_fig.hPlot_traces_thresh, ...
        S_fig.hPlot_traces_spk, S_fig.hPlot_traces_spk1);
        title_(S_fig.hAx_traces, sprintf('Scale: %0.1f uV', S_fig.maxAmp));
        set(hFig, 'UserData', S_fig);

        case {'leftarrow', 'rightarrow', 'home', 'end'} %navigation
        switch lower(event.Key)
            case 'leftarrow'
            nlim_bin = S_fig.nlim_bin - (S_fig.nLoad_bin) * factor; %no overlap
            if nlim_bin(1)<1
                msgbox_('Beginning of file', 1);
                nlim_bin = [1, S_fig.nLoad_bin];
            end
            case 'rightarrow'
            nlim_bin = S_fig.nlim_bin + (S_fig.nLoad_bin + 1) * factor; %no overlap
            if nlim_bin(2) > S_fig.nSamples_bin
                msgbox_('End of file', 1);
                nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
            end
            case 'home', nlim_bin = [1, S_fig.nLoad_bin]; %begging of file
            case 'end', nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin; %end of file
        end %switch
        S_fig.nlim_bin = nlim_bin;
        set(hFig, 'UserData', S_fig);
        Fig_preview_plot_(P);

        case 'f' %apply filter
        set(hFig, 'UserData', setfield(S_fig, 'fFilter', ~S_fig.fFilter));
        Fig_preview_plot_(P, 1);

        case 't' %toggle spike threshold
        set(hFig, 'UserData', setfield(S_fig, 'fThresh_spk', ~S_fig.fThresh_spk));
        Fig_preview_plot_(P, 1);

        case 's' %show/hide spikes
        set(hFig, 'UserData', setfield(S_fig, 'fShow_spk', ~S_fig.fShow_spk));
        Fig_preview_plot_(P, 1);

        case 'g' %grid toggle on/off
        S_fig.fGrid = ~S_fig.fGrid;
        grid_([S_fig.hAx_traces, S_fig.hAx_mean, S_fig.hAx_psd, S_fig.hAx_sites], S_fig.fGrid);
        set(hFig, 'UserData', S_fig);
        menu_label_('menu_preview_view_grid', ifeq_(S_fig.fGrid, 'Hide [G]rid', 'Show [G]rid'));

        case 'h', msgbox_(S_fig.csHelp, 1); %help

        case 'r' %reset view
        P = get0_('P');
        axis_(S_fig.hAx_traces, [S_fig.nlim_bin / P.sampleRateHz, S_fig.siteLim+[-1,1]]);

        case 'e' % export to workspace
        Fig_preview_export_(hFig);

    end %switch
end % function
