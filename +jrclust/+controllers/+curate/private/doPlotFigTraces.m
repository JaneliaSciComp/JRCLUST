function tracesFilt = doPlotFigTraces(hFigTraces, hCfg, tracesRaw, resetAxis, hClust)
    % fAxis_reset: reset the axis limit
    % [usage]
    % doPlotFigTraces()
    % doPlotFigTraces(fAxis_reset)
    % 2017/06/22 James Jun
    % 6/22 JJJ: added seperator lines, fixed the reset view and spike view

    if nargin < 4
        resetAxis = 0;
    end
    if nargin < 5
        hClust = [];
    end

    fWait = msgbox_('Plotting...', 0, 1);
    fShuttleOrder = 1; % shuffle cluster color

    S0 = get0_();

    hFigTraces.wait(true);

    sampleRate = hCfg.sampleRate / hCfg.nSkip_show;
    viSamples1 = 1:hCfg.nSkip_show:size(tracesRaw, 1);
    evtWindowSamp = round(hCfg.evtWindowSamp / hCfg.nSkip_show); %show 2x of range

    if strcmpi(hFigTraces.figData.filter, 'on')
        hCfg.sampleRate = sampleRate;
        hCfg.useGPU = false;
        hCfg.filterType = hCfg.dispFilter;

        if hCfg.fft_thresh > 0
            tracesRaw = jrclust.filters.fftClean(tracesRaw, hCfg);
        end

        tracesFilt = jrclust.utils.bit2uV(jrclust.filters.filtCAR(tracesRaw(viSamples1, :), hCfg), hCfg);
        vcFilter_show = hCfg.filterType;
    else
        tracesFilt = jrclust.utils.meanSubtract(single(tracesRaw(viSamples1, :))) * hCfg.bitScaling;
        vcFilter_show = 'off';
    end

    viSites = 1:numel(hCfg.siteMap);

    nTime_traces = get_(hCfg, 'nTime_traces');
    if isempty(nTime_traces) || nTime_traces==1
        vrTime_bin = ((hFigTraces.figData.nlim_bin(1):hCfg.nSkip_show:hFigTraces.figData.nlim_bin(end))-1) / hCfg.sampleRate;
        vcXLabel = 'Time (s)';
    else
        vrTime_bin = (0:(size(tracesFilt,1)-1)) / (hCfg.sampleRate / hCfg.nSkip_show) + (hFigTraces.figData.nlim_bin(1)-1) / hCfg.sampleRate;
        [cvn_lim_bin, viRange_bin, viEdges] = sample_skip_(hFigTraces.figData.nlim_bin, hFigTraces.figData.nSamples_bin, nTime_traces);
        tlim_show = (cellfun(@(x)x(1), cvn_lim_bin([1,end]))) / hCfg.sampleRate;
        vcXLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', nTime_traces, tlim_show, diff(hCfg.tlim));
        mrX_edges = vrTime_bin(repmat(viEdges(:)', [3,1]));
        mrY_edges = repmat([0;numel(hCfg.siteMap)+1;nan],1,numel(viEdges));
        set(hFigTraces.figData.hPlot_edges, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
        csTime_bin = cellfun(@(x)sprintf('%0.1f', x(1)/hCfg.sampleRate), cvn_lim_bin, 'UniformOutput', 0);
        set(hFigTraces.figData.hAx, {'XTick', 'XTickLabel'}, {vrTime_bin(viEdges), csTime_bin});
    end
    multiplot(hFigTraces.figData.hPlot, hFigTraces.figData.maxAmp, vrTime_bin, tracesFilt, viSites);
    % axis(S_fig.hAx, [vrTime_bin(1), vrTime_bin(end), viSites(1)-1, viSites(end)+1]);
    grid(hFigTraces.figData.hAx, hFigTraces.figData.vcGrid);
    set(hFigTraces.figData.hAx, 'YTick', viSites);
    title_(hFigTraces.figData.hAx, sprintf(hFigTraces.figData.vcTitle, hFigTraces.figData.maxAmp));
    xlabel(hFigTraces.figData.hAx, vcXLabel);
    ylabel(hFigTraces.figData.hAx, 'Site #');
    set(hFigTraces.figData.hPlot, 'Visible', hFigTraces.figData.vcTraces);

    % Delete spikes from other threads
    S_fig_ = get(hFigTraces, 'UserData');
    if isfield(S_fig_, 'chSpk'), delete_multi_(S_fig_.chSpk); end
    if isfield(hFigTraces.figData, 'chSpk'), delete_multi_(hFigTraces.figData.chSpk); end

    % plot spikes
    if strcmpi(hFigTraces.figData.vcSpikes, 'on') && isfield(S0, 'viTime_spk')
        viTime_spk = S0.viTime_spk - int32(S0.viT_offset_file(S0.iFile_show));
        if nTime_traces > 1
            viSpk1 = find(in_range_(viTime_spk, cvn_lim_bin));
            [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
            viTime_spk1 = round(reverse_lookup_(viTime_spk1, viRange_bin) / hCfg.nSkip_show);
        else
            viSpk1 = find(viTime_spk >= hFigTraces.figData.nlim_bin(1) & viTime_spk < hFigTraces.figData.nlim_bin(end));
            [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
            viTime_spk1 = round((viTime_spk1 - hFigTraces.figData.nlim_bin(1) + 1) / hCfg.nSkip_show); %time offset
        end
        t_start1 = single(hFigTraces.figData.nlim_bin(1) - 1) / hCfg.sampleRate;
        viSite_spk1 = single(viSite_spk1);
        % check if clustered
        if isempty(hClust)
            nSites = size(tracesFilt,2);
            chSpk = cell(nSites, 1);
            for iSite=1:nSites %deal with subsample factor
                viSpk11 = find(viSite_spk1 == iSite);
                if isempty(viSpk11), continue; end
                viTime_spk11 = viTime_spk1(viSpk11);
                [mrY11, mrX11] = vr2mr3_(tracesFilt(:,iSite), viTime_spk11, evtWindowSamp); %display purpose x2
                %             vr2mr_spk_(tracesFilt(:,iSite), viTime_spk11, P);
                mrT11 = single(mrX11-1) / sampleRate + t_start1;
                chSpk{iSite} = line(nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5, 'Parent', hFigTraces.figData.hAx);
                multiplot(chSpk{iSite}, hFigTraces.figData.maxAmp, mrT11, mrY11, iSite);
            end
        else % different color for each clu
            viClu_spk1 = hClust.viClu(viSpk1);
            mrColor_clu = [jet(hClust.nClu); 0 0 0];
            vrLineWidth_clu = (mod((1:hClust.nClu)-1, 3)+1)'/2 + .5;  %(randi(3, S_clu.nClu, 1)+1)/2;
            if fShuttleOrder
                mrColor_clu = shuffle_static_(mrColor_clu, 1);
                vrLineWidth_clu = shuffle_static_(vrLineWidth_clu, 1);
            end
            nSpk1 = numel(viTime_spk1);
            chSpk = cell(nSpk1, 1);
            for iSpk1 = 1:nSpk1
                iTime_spk11 = viTime_spk1(iSpk1);
                iSite11 = viSite_spk1(iSpk1);
                [mrY11, mrX11] = vr2mr3_(tracesFilt(:,iSite11), iTime_spk11, evtWindowSamp); %display purpose x2
                mrT11 = double(mrX11-1) / sampleRate + t_start1;
                iClu11 = viClu_spk1(iSpk1);
                if iClu11<=0, continue; end
                %                 vrColor1 = [0 0 0];
                %                 lineWidth1 = .5;
                %             else
                vrColor1 = mrColor_clu(iClu11,:);
                lineWidth1 = vrLineWidth_clu(iClu11);
                %             end
                chSpk{iSpk1} = line(nan, nan, 'Color', vrColor1, 'LineWidth', lineWidth1, 'Parent', hFigTraces.figData.hAx);
                multiplot(chSpk{iSpk1}, hFigTraces.figData.maxAmp, mrT11, mrY11, iSite11);
            end
        end
        hFigTraces.figData.chSpk = chSpk;
    else
        % delete spikes
        hFigTraces.figData.chSpk = [];
    end
    if resetAxis, fig_traces_reset_(hFigTraces.figData); end
    set(hFigTraces, 'UserData', hFigTraces.figData, 'Name', sprintf('%s: filter: %s', hCfg.vcFile_prm, (vcFilter_show)));
    figure_wait_(0, hFigTraces);
    close_(fWait);
end %func
