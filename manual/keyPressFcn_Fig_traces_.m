%--------------------------------------------------------------------------
function keyPressFcn_Fig_traces_(hFig, event)
    % 2017/6/22 James Jun: Added nTime_traces multiview

    global mnWav1 mrWav1 mnWav
    S0 = get(0, 'UserData');
    P = S0.P;
    S_fig = get(hFig, 'UserData');
    factor = 1 + 3 * keyModifier(event, 'shift');
    nSites = numel(P.chanMap);

    switch lower(event.Key)
        case 'h', msgbox_(S_fig.csHelp, 1);

        case {'uparrow', 'downarrow'}
        if isfield(S_fig, 'chSpk')
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot, S_fig.chSpk);
        else
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot);
        end
        title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));
        set(hFig, 'UserData', S_fig);

        case {'leftarrow', 'rightarrow', 'j', 'home', 'end'}
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
            case 'home' %beginning of file
            nlim_bin = [1, S_fig.nLoad_bin];
            case 'end' %end of file
            nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
            case 'j'
            vcAns = inputdlg_('Go to time (s)', 'Jump to time', 1, {'0'});
            if isempty(vcAns), return; end
            try
                nlim_bin = round(str2double(vcAns)*P.sRateHz) + [1, S_fig.nLoad_bin];
            catch
                return;
            end
        end %switch
        nTime_traces = get_(P, 'nTime_traces');
        [cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, S_fig.nSamples_bin, nTime_traces);
        if P.fTranspose_bin
            fseek_(S_fig.fid_bin, nlim_bin(1), P);
            if nTime_traces > 1
                mnWav1 = load_bin_multi_(S_fig.fid_bin, cvn_lim_bin, P)';
            else
                mnWav1 = load_bin_(S_fig.fid_bin, P.vcDataType, [P.nChans, S_fig.nLoad_bin])';
            end
        else
            mnWav1 = mnWav(viRange_bin, :);
        end
        mnWav1 = uint2int_(mnWav1);
        S_fig.nlim_bin = nlim_bin;
        set_fig_(hFig, S_fig);
        Fig_traces_plot_(1); %redraw

        case 'f' %apply filter
        S_fig.vcFilter = str_toggle_(S_fig.vcFilter, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();

        case 'g' %grid toggle on/off
        S_fig.vcGrid = str_toggle_(S_fig.vcGrid, 'on', 'off');
        grid(S_fig.hAx, S_fig.vcGrid);
        set(hFig, 'UserData', S_fig);

        case 'r' %reset view
        fig_traces_reset_(S_fig);

        case 'e' %export current view
        assignWorkspace_(mnWav1, mrWav1);
        disp('mnWav1: raw traces, mrWav1: filtered traces');

        case 'p' %power spectrum
        iSite_show = inputdlg_num_(sprintf('Site# to show (1-%d, 0 for all)', nSites), 'Site#', 0);
        if isnan(iSite_show), return; end
        hFig = createFigure('FigPsd', [.5 0 .5 1], P.paramFile, 1, 1); %show to the right
        % ask user which channels to plot
        if iSite_show>0
            mrWav2 = mrWav1(:, iSite_show);
        else
            mrWav2 = mrWav1;
        end
        plotMedPower_(mrWav2, 'sRateHz', P.sRateHz/P.nSkip_show, 'viChanExcl', P.viSiteZero);

        case 's' %show/hide spikes
        S_fig.vcSpikes = str_toggle_(S_fig.vcSpikes, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();

        case 't' %show/hide traces
        S_fig.vcTraces = str_toggle_(S_fig.vcTraces, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();

        case 'c' %channel query
        msgbox_('Draw a rectangle', 1);
        hRect = imrect_();
        if isempty(hRect), return ;end
        vrPos_rect = getPosition(hRect);
        S_plot = get(S_fig.hPlot, 'UserData');
        vrX = get(S_fig.hPlot, 'XData');
        vrY = get(S_fig.hPlot, 'YData');
        viIndex = find(vrX >= vrPos_rect(1) & vrX <= sum(vrPos_rect([1,3])) & vrY >= vrPos_rect(2) & vrY <= sum(vrPos_rect([2,4])));
        if isempty(viIndex), deleteMany(hRect); return; end
        index_plot = round(median(viIndex));
        [time1, iSite] = ind2sub(size(mrWav1), index_plot);
        mrX = reshape(vrX, S_plot.dimm);
        mrY = reshape(vrY, S_plot.dimm);
        hold(S_fig.hAx, 'on');
        hPoint = plot(vrX(index_plot), vrY(index_plot), 'r*');
        hLine = plot(S_fig.hAx, mrX(:,iSite), mrY(:,iSite), 'r-');
        hold(S_fig.hAx, 'off');
        iChan = P.chanMap(iSite);
        msgbox_(sprintf('Site: %d/ Chan: %d', iSite, iChan), 1);
        deleteMany(hRect, hLine, hPoint);
    end %return if S_fig didn't change
end %func
