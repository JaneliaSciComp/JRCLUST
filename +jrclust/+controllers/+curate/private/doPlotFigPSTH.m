function [hFigTrial1, hFigTrial2] = doPlotFigPSTH(hClust, hFigTrial1, hFigTrial2, selected)
    %doPlotFigPSTH()
    %   plot if window open using curretnly selected clusters
    %doPlotFigPSTH(hCfg, iClu, S_clu)
    %   Open window and plot specific clusters and S_clu
    hCfg = hClust.hCfg;

    % begin TW block
    if isempty(hCfg.trialFile)
        if exist(jrclust.utils.subsExt(hCfg.configFile, '.starts.mat'), 'file')
            hCfg.trialFile = jrclust.utils.subsExt(hCfg.configFile, '.starts.mat');
        elseif exist(jrclust.utils.subsExt(hCfg.configFile, '.mat'), 'file')
            hCfg.trialFile = jrclust.utils.subsExt(hCfg.configFile, '.mat');
        else
            jrclust.utils.qMsgBox('''trialFile'' not set. Reload .prm file after setting (under "File menu")');
            return;
        end
    end
    % end TW block

    % import trial times
    trialTimes = loadTrialFile(hCfg.trialFile);
    if ~iscell(trialTimes)
        trialTimes = {trialTimes};
    end

    if isempty(trialTimes) % failed to load
        jrclust.utils.qMsgBox('Trial file does not exist', 0, 1);
        return;
    end

    nStims = numel(trialTimes);

    %% plot primary/secondary figures
    axOffset = 0.08;
    axLen = 1/nStims;

    if ~jrclust.utils.isvalid(hFigTrial1) || ~hFigTrial1.isReady
        hFigTrial1 = jrclust.views.Figure('FigTrial1', [.5  .5 .5 .5], hCfg.trialFile, 0, 0);
%         [vhAx1, vhAx2] = deal(nan(nStims, 1));
        for iStim = 1:nStims
            iOffset = axOffset + (iStim-1) * axLen;
            hFigTrial1.addAxes(sprintf('stim%d1', iStim), 'Position', [.08 iOffset .9 axLen*.68]);
            hFigTrial1.addAxes(sprintf('stim%d2', iStim), 'Position', [.08 iOffset + axLen*.68 .9 axLen*.2]);
            % vhAx1(iStim) = axes('Parent', hFigTrial1, 'Position', [.08 iOffset .9 axLen*.68]);
            % vhAx2(iStim) = axes('Parent', hFigTrial1, 'Position', [.08 iOffset + axLen*.68 .9 axLen*.2]);
        end
        hFigTrial1.figData.color = 'k';
%         vcColor = 'k';
%         set(hFigTrial1, 'UserData', makeStruct_(vhAx1, vhAx2, vcColor));
    end
    plot_figure_psth_(hFigTrial1, selected(1), trialTimes, hClust, hCfg);

    if ~jrclust.utils.isvalid(hFigTrial2) || ~hFigTrial2.isReady
        hFigTrial2 = jrclust.views.Figure('FigTrial2', [.5  0 .5 .5], hCfg.trialFile, 0, 0);
        hFigTrial2.figApply(@set, 'Visible', 'off');
%         [vhAx1, vhAx2] = deal(nan(nStims, 1));
        for iStim = 1:nStims
            iOffset = axOffset + (iStim-1) * axLen;
            hFigTrial2.addAxes(sprintf('stim%d1', iStim), 'Position', [.08 iOffset .9 axLen*.68]);
            hFigTrial2.addAxes(sprintf('stim%d2', iStim), 'Position', [.08 iOffset + axLen*.68 .9 axLen*.2]);
            % vhAx1(iStim) = axes('Parent', hFigTrial1, 'Position', [.08 iOffset .9 axLen*.68]);
            % vhAx2(iStim) = axes('Parent', hFigTrial1, 'Position', [.08 iOffset + axLen*.68 .9 axLen*.2]);
        end
        hFigTrial2.figData.color = 'r';
%         vcColor = 'r';
%         set(hFigTrial2, 'UserData', makeStruct_(vhAx1, vhAx2, vcColor));
    end
    if numel(selected) == 2
        hFigTrial2.figApply(@set, 'Visible', 'on');
        plot_figure_psth_(hFigTrial2, selected(2), trialTimes, hClust, hCfg);
    else
        hFigTrial2.figApply(@set, 'Visible', 'off');
    end
end

%% LOCAL FUNCTIONS
function trialTimes = loadTrialFile(trialFile)
    %LOADTRIALFILE Import trial times (in seconds)
    trialTimes = [];

    try
        [~, ~, ext] = fileparts(trialFile);

        if strcmpi(ext, '.mat')
            trialData = load(trialFile);
            fieldNames = fieldnames(trialData);

            trialTimes = trialData.(fieldNames{1});
            if isstruct(trialTimes)
                trialTimes = trialTimes.times;
            end
        elseif strcmpi(ext, '.csv')
            trialTimes = csvread(trialFile);
            if isrow(trialTimes)
                trialTimes = trialTimes(:);
            end
        end
    catch ME
        warning(ME.identifier, 'Could not load trialFile %s: %s', trialFile, ME.message);
    end
end

function plot_figure_psth_(hFigTrial, iCluster, trialTimes, hClust, hCfg)
%     [vhAx1, vhAx2] = deal(S_fig.vhAx1, S_fig.vhAx2, S_fig.vcColor);
    hAxes = keys(hFigTrial.hAxes);
    nStims = numel(hAxes)/2;

    for iStim = 1:nStims
        axKey1 = sprintf('stim%d1', iStim); hAx1 = hFigTrial.hAxes(axKey1);
        axKey2 = sprintf('stim%d2', iStim); hAx2 = hFigTrial.hAxes(axKey2);

        cla(hAx1);
        cla(hAx2);
        iTrialTimes = trialTimes{iStim}; %(:,1);
        nTrials = numel(iTrialTimes);
        clusterTimes = hClust.spikeTimes(hClust.spikesByCluster{iCluster});
        plot_raster_clu_(clusterTimes, iTrialTimes, hCfg, hAx1);
        plot_psth_clu_(clusterTimes, iTrialTimes, hCfg, hAx2, hFigTrial.figData.color);
        hFigTrial.axApply(axKey2, @title, sprintf('Cluster %d; %d trials', iCluster, nTrials));
    end
    %     offset = offset + nTrials;
    if nStims > 1
        arrayfun(@(i) hFigTrial.axApply(sprintf('stim%d1', i), @set, 'XTickLabel', {}), 2:nStims);
        arrayfun(@(i) hFigTrial.axApply(sprintf('stim%d1', i), @xlabel, ''), 2:nStims);
        %set(vhAx1(2:end),'xticklabel',{});
%         for ax = vhAx1(2:end)
%             xlabel(ax, '')
%         end
    end
end

function plot_raster_clu_(clusterTimes, trialTimes, hCfg, hAx)
    trialLength = diff(hCfg.tlim_psth); % seconds
    nTrials = numel(trialTimes);
    spikeTimes = cell(nTrials, 1);
    t0 = -hCfg.tlim_psth(1);
    for iTrial = 1:nTrials
        rTime_trial1 = trialTimes(iTrial);
        vrTime_lim1 = rTime_trial1 + hCfg.tlim_psth;
        vrTime_clu1 = double(clusterTimes) / hCfg.sampleRate;
        vrTime_clu1 = vrTime_clu1(vrTime_clu1>=vrTime_lim1(1) & vrTime_clu1<vrTime_lim1(2));
        vrTime_clu1 = (vrTime_clu1 - rTime_trial1 + t0) / trialLength;
        spikeTimes{iTrial} = vrTime_clu1';
    end

    % Plot
    plotSpikeRaster(spikeTimes,'PlotType','vertline','RelSpikeStartTime',0,'XLimForCell',[0 1], ...
        'LineFormat', struct('LineWidth', 1.5), 'hAx', hAx);
    ylabel(hAx, 'Trial #')
    % title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
    vrXTickLabel = hCfg.tlim_psth(1):(hCfg.xtick_psth):hCfg.tlim_psth(2);
    vrXTick = linspace(0,1,numel(vrXTickLabel));
    set(hAx, {'XTick', 'XTickLabel'}, {vrXTick, vrXTickLabel});
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [t0,t0]/trialLength, get(hAx,'YLim'), 'r-');
    xlabel(hAx, 'Time (s)');
end

function plot_psth_clu_(clusterTimes, trialTimes, hCfg, hAx, vcColor)
    tbin = hCfg.tbin_psth;
    nbin = round(tbin * hCfg.sampleRate);
    nlim = round(hCfg.tlim_psth/tbin);
    viTime_Trial = round(trialTimes / tbin);

    vlTime1 = zeros(0);
    vlTime1(ceil(double(clusterTimes)/nbin)) = 1;
    mr1 = vr2mr2_(double(vlTime1), viTime_Trial, nlim);
    vnRate = mean(mr1,2) / tbin;
    vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;
    bar(hAx, vrTimePlot, vnRate, 1, 'EdgeColor', 'none', 'FaceColor', vcColor);
    vrXTick = hCfg.tlim_psth(1):(hCfg.xtick_psth):hCfg.tlim_psth(2);
    set(hAx, 'XTick', vrXTick, 'XTickLabel', []);
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [0 0], get(hAx, 'YLim'), 'r-');
    ylabel(hAx, 'Rate (Hz)');
    xlim(hAx, hCfg.tlim_psth);
end
