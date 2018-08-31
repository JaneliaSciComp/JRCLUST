%--------------------------------------------------------------------------
% plot data
function Fig_traces_plot_(fAxis_reset)
    % fAxis_reset: reset the axis limit
    % [usage]
    % Fig_traces_plot_()
    % Fig_traces_plot_(fAxis_reset)
    % 2017/06/22 James Jun
    % 6/22 JJJ: added seperator lines, fixed the reset view and spike view

    global mnWav1 mrWav1 % current timeslice to plot
    if nargin<1, fAxis_reset = 0; end
    fWait = msgbox_('Plotting...',0,1);
    fShuttleOrder = 1; %shuffle cluster color
    [S0, P, S_clu] = get0_();
    [hFig, S_fig] = getCachedFig('Fig_traces');
    figure_wait_(1, hFig); drawnow;
    sampleRateHz = P.sampleRateHz / P.nSkip_show;
    viSamples1 = 1:P.nSkip_show:size(mnWav1,1);
    spkLim = round(P.spkLim / P.nSkip_show); %show 2x of range
    P.vcFilter = get_filter_(P);
    if strcmpi(S_fig.vcFilter, 'on')
        P1=P; P1.sampleRateHz = sampleRateHz; P1.useGPU = 0;
        P1.vcFilter = getOr(P, 'vcFilter_show', P.vcFilter);
        if P.fft_thresh>0, mnWav1 = fft_clean_(mnWav1, P); end
        mrWav1 = bit2uV_(filt_car_(mnWav1(viSamples1, P.chanMap), P1), P1);
        vcFilter_show = P1.vcFilter;
    else
        mrWav1 = meanSubtract(single(mnWav1(viSamples1, P.chanMap))) * P.uV_per_bit;
        vcFilter_show = 'off';
    end
    viSites = 1:numel(P.chanMap);
    % mrWav1 = meanSubtract(single(mnWav1(:, P.chanMap))) * P.uV_per_bit;
    % hide bad channels
    nTime_traces = get_(P, 'nTime_traces');
    if isempty(nTime_traces) || nTime_traces==1
        vrTime_bin = ((S_fig.nlim_bin(1):P.nSkip_show:S_fig.nlim_bin(end))-1) / P.sampleRateHz;
        vcXLabel = 'Time (s)';
    else
        vrTime_bin = (0:(size(mrWav1,1)-1)) / (P.sampleRateHz / P.nSkip_show) + (S_fig.nlim_bin(1)-1) / P.sampleRateHz;
        [cvn_lim_bin, viRange_bin, viEdges] = sample_skip_(S_fig.nlim_bin, S_fig.nSamples_bin, nTime_traces);
        tlim_show = (cellfun(@(x)x(1), cvn_lim_bin([1,end]))) / P.sampleRateHz;
        vcXLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', nTime_traces, tlim_show, diff(P.tlim));
        mrX_edges = vrTime_bin(repmat(viEdges(:)', [3,1]));
        mrY_edges = repmat([0;numel(P.chanMap)+1;nan],1,numel(viEdges));
        set(S_fig.hPlot_edges, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
        csTime_bin = cellfun(@(x)sprintf('%0.1f', x(1)/P.sampleRateHz), cvn_lim_bin, 'UniformOutput', 0);
        set(S_fig.hAx, {'XTick', 'XTickLabel'}, {vrTime_bin(viEdges), csTime_bin});
    end
    multiplot(S_fig.hPlot, S_fig.maxAmp, vrTime_bin, mrWav1, viSites);
    % axis(S_fig.hAx, [vrTime_bin(1), vrTime_bin(end), viSites(1)-1, viSites(end)+1]);
    grid(S_fig.hAx, S_fig.vcGrid);
    set(S_fig.hAx, 'YTick', viSites);
    title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));
    xlabel(S_fig.hAx, vcXLabel);
    ylabel(S_fig.hAx, 'Site #');
    set(S_fig.hPlot, 'Visible', S_fig.vcTraces);

    % Delete spikes from other threads
    S_fig_ = get(hFig, 'UserData');
    if isfield(S_fig_, 'chSpk'), deleteMany(S_fig_.chSpk); end
    if isfield(S_fig, 'chSpk'), deleteMany(S_fig.chSpk); end

    % plot spikes
    if strcmpi(S_fig.vcSpikes, 'on') && isfield(S0, 'spikeTimes')
        spikeTimes = S0.spikeTimes - int32(S0.fileSampleOffsets(S0.iFile_show));
        if nTime_traces > 1
            viSpk1 = find(in_range_(spikeTimes, cvn_lim_bin));
            [spikeSites1, spikeTimes1] = multifun_(@(vr)vr(viSpk1), S0.spikeSites, spikeTimes);
            spikeTimes1 = round(reverse_lookup_(spikeTimes1, viRange_bin) / P.nSkip_show);
        else
            viSpk1 = find(spikeTimes >= S_fig.nlim_bin(1) & spikeTimes < S_fig.nlim_bin(end));
            [spikeSites1, spikeTimes1] = multifun_(@(vr)vr(viSpk1), S0.spikeSites, spikeTimes);
            spikeTimes1 = round((spikeTimes1 - S_fig.nlim_bin(1) + 1) / P.nSkip_show); %time offset
        end
        t_start1 = single(S_fig.nlim_bin(1) - 1) / P.sampleRateHz;
        spikeSites1 = single(spikeSites1);
        % check if clustered
        if isempty(S_clu)
            nSites = size(mrWav1,2);
            chSpk = cell(nSites, 1);
            for iSite=1:nSites %deal with subsample factor
                viSpk11 = find(spikeSites1 == iSite);
                if isempty(viSpk11), continue; end
                spikeTimes11 = spikeTimes1(viSpk11);
                [mrY11, mrX11] = vr2mr3_(mrWav1(:,iSite), spikeTimes11, spkLim); %display purpose x2
                %             vr2mr_spk_(mrWav1(:,iSite), spikeTimes11, P);
                mrT11 = single(mrX11-1) / sampleRateHz + t_start1;
                chSpk{iSite} = line(nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5, 'Parent', S_fig.hAx);
                multiplot(chSpk{iSite}, S_fig.maxAmp, mrT11, mrY11, iSite);
            end
        else % different color for each clu
            viClu_spk1 = S_clu.spikeClusters(viSpk1);
            mrColor_clu = [jet(S_clu.nClusters); 0 0 0];
            vrLineWidth_clu = (mod((1:S_clu.nClusters)-1, 3)+1)'/2 + .5;  %(randi(3, S_clu.nClusters, 1)+1)/2;
            if fShuttleOrder
                mrColor_clu = shuffle_static_(mrColor_clu, 1);
                vrLineWidth_clu = shuffle_static_(vrLineWidth_clu, 1);
            end
            nSpk1 = numel(spikeTimes1);
            chSpk = cell(nSpk1, 1);
            for iSpk1 = 1:nSpk1
                iTime_spk11 = spikeTimes1(iSpk1);
                iSite11 = spikeSites1(iSpk1);
                [mrY11, mrX11] = vr2mr3_(mrWav1(:,iSite11), iTime_spk11, spkLim); %display purpose x2
                mrT11 = double(mrX11-1) / sampleRateHz + t_start1;
                iClu11 = viClu_spk1(iSpk1);
                if iClu11<=0, continue; end
                %                 vrColor1 = [0 0 0];
                %                 lineWidth1 = .5;
                %             else
                vrColor1 = mrColor_clu(iClu11,:);
                lineWidth1 = vrLineWidth_clu(iClu11);
                %             end
                chSpk{iSpk1} = line(nan, nan, 'Color', vrColor1, 'LineWidth', lineWidth1, 'Parent', S_fig.hAx);
                multiplot(chSpk{iSpk1}, S_fig.maxAmp, mrT11, mrY11, iSite11);
            end
        end
        S_fig.chSpk = chSpk;
    else
        % delete spikes
        S_fig.chSpk = [];
    end
    if fAxis_reset, fig_traces_reset_(S_fig); end
    set(hFig, 'UserData', S_fig, 'Name', sprintf('%s: filter: %s', P.paramFile, (vcFilter_show)));
    figure_wait_(0, hFig);
    tryClose(fWait);
end % function
