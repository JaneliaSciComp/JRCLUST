%--------------------------------------------------------------------------
function plot_raster_clu_(viTime_clu, vrTime_trial, P, hAx)
    if nargin<4, hAx=gca; end

    trialLength = diff(P.tlim_psth); % seconds
    nTrials = numel(vrTime_trial);
    spikeTimes = cell(nTrials, 1);
    t0 = -P.tlim_psth(1);
    for iTrial = 1:nTrials
        rTime_trial1 = vrTime_trial(iTrial);
        vrTime_lim1 = rTime_trial1 + P.tlim_psth;
        vrTime_clu1 = double(viTime_clu) / P.sampleRateHz;
        vrTime_clu1 = vrTime_clu1(vrTime_clu1>=vrTime_lim1(1) & vrTime_clu1<vrTime_lim1(2));
        vrTime_clu1 = (vrTime_clu1 - rTime_trial1 + t0) / trialLength;
        spikeTimes{iTrial} = vrTime_clu1';
    end

    % Plot
    % hAx_pre = axes_(hAx);
    plotSpikeRaster(spikeTimes,'PlotType','vertline','RelSpikeStartTime',0,'XLimForCell',[0 1], ...
    'LineFormat', struct('LineWidth', 1.5), 'hAx', hAx);
    % axes_(hAx_pre);
    ylabel(hAx, 'Trial #')
    % title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
    vrXTickLabel = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
    vrXTick = linspace(0,1,numel(vrXTickLabel));
    set(hAx, {'XTick', 'XTickLabel'}, {vrXTick, vrXTickLabel});
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [t0,t0]/trialLength, get(hAx,'YLim'), 'r-');
    xlabel(hAx, 'Time (s)');
end % function
