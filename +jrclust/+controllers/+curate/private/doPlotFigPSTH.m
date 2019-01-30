function [hFigTrial1, hFigTrial2] = doPlotFigPSTH(hClust, hFigTrial1, hFigTrial2, selected)
    %DOPLOTFIGPSTH Plot PSTH figures
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
        for iStim = 1:nStims
            iOffset = axOffset + (iStim-1) * axLen;
            hFigTrial1.addAxes(sprintf('stim%d1', iStim), 'Position', [axOffset iOffset .9 axLen*.68]);
            hFigTrial1.addAxes(sprintf('stim%d2', iStim), 'Position', [axOffset iOffset + axLen*.68 .9 axLen*.2]);
        end

        hFigTrial1.figData.color = 'k';
    end
    plot_figure_psth_(hFigTrial1, selected(1), trialTimes, hClust, hCfg);

    if ~jrclust.utils.isvalid(hFigTrial2) || ~hFigTrial2.isReady
        hFigTrial2 = jrclust.views.Figure('FigTrial2', [.5  0 .5 .5], hCfg.trialFile, 0, 0);
        hFigTrial2.figApply(@set, 'Visible', 'off');
        for iStim = 1:nStims
            iOffset = axOffset + (iStim-1) * axLen;
            hFigTrial2.addAxes(sprintf('stim%d1', iStim), 'Position', [axOffset iOffset .9 axLen*.68]);
            hFigTrial2.addAxes(sprintf('stim%d2', iStim), 'Position', [axOffset iOffset + axLen*.68 .9 axLen*.2]);
        end

        hFigTrial2.figData.color = 'r';
    end

    % show this plot iff we have a second selected cluster
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
        warning('Could not load trialFile %s: %s', trialFile, ME.message);
    end
end

function plot_figure_psth_(hFigTrial, iCluster, trialTimes, hClust, hCfg)
%     [vhAx1, vhAx2] = deal(S_fig.vhAx1, S_fig.vhAx2, S_fig.vcColor);
    hAxes = keys(hFigTrial.hAxes);
    nStims = numel(hAxes)/2;

    for iStim = 1:nStims
        axKey1 = sprintf('stim%d1', iStim); hAx1 = hFigTrial.hAxes(axKey1); cla(hAx1);
        axKey2 = sprintf('stim%d2', iStim); hAx2 = hFigTrial.hAxes(axKey2); cla(hAx2);

        iTrialTimes = trialTimes{iStim}; %(:,1);
        nTrials = numel(iTrialTimes);
        clusterTimes = hClust.spikeTimes(hClust.spikesByCluster{iCluster});
        plot_raster_clu_(clusterTimes, iTrialTimes, hCfg, hAx1);
        plot_psth_clu_(clusterTimes, iTrialTimes, hCfg, hAx2, hFigTrial.figData.color);
        hFigTrial.axApply(axKey2, @title, sprintf('Cluster %d; %d trials', iCluster, nTrials));
    end

    if nStims > 1
        arrayfun(@(i) hFigTrial.axApply(sprintf('stim%d1', i), @set, 'XTickLabel', {}), 2:nStims);
        arrayfun(@(i) hFigTrial.axApply(sprintf('stim%d1', i), @xlabel, ''), 2:nStims);
    end
end

function plot_raster_clu_(clusterTimes, trialTimes, hCfg, hAx)
    trialLength = diff(hCfg.psthTimeLimits); % seconds
    nTrials = numel(trialTimes);
    spikeTimes = cell(nTrials, 1);
    t0 = -hCfg.psthTimeLimits(1);
    for iTrial = 1:nTrials
        rTime_trial1 = trialTimes(iTrial);
        vrTime_lim1 = rTime_trial1 + hCfg.psthTimeLimits;
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
    vrXTickLabel = hCfg.psthTimeLimits(1):(hCfg.psthXTick):hCfg.psthTimeLimits(2);
    vrXTick = linspace(0,1,numel(vrXTickLabel));
    set(hAx, {'XTick', 'XTickLabel'}, {vrXTick, vrXTickLabel});
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [t0,t0]/trialLength, get(hAx,'YLim'), 'r-');
    xlabel(hAx, 'Time (s)');
end

function plot_psth_clu_(clusterTimes, trialTimes, hCfg, hAx, vcColor)
    tbin = hCfg.psthTimeBin;
    nbin = round(tbin * hCfg.sampleRate);
    nlim = round(hCfg.psthTimeLimits/tbin);
    viTime_Trial = round(trialTimes / tbin);

    vlTime1 = zeros(0);
    vlTime1(ceil(double(clusterTimes)/nbin)) = 1;
    mr1 = vr2mr2_(double(vlTime1), viTime_Trial, nlim);
    vnRate = mean(mr1,2) / tbin;
    vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;
    bar(hAx, vrTimePlot, vnRate, 1, 'EdgeColor', 'none', 'FaceColor', vcColor);
    vrXTick = hCfg.psthTimeLimits(1):(hCfg.psthXTick):hCfg.psthTimeLimits(2);
    set(hAx, 'XTick', vrXTick, 'XTickLabel', []);
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [0 0], get(hAx, 'YLim'), 'r-');
    ylabel(hAx, 'Rate (Hz)');
    xlim(hAx, hCfg.psthTimeLimits);
end

function mr = vr2mr2_(vr, viRow, spkLim, viCol)
    if nargin<4, viCol = []; end
    % JJJ 2015 Dec 24
    % vr2mr2: quick version and doesn't kill index out of range
    % assumes vi is within range and tolerates spkLim part of being outside
    % works for any datatype

    % prepare indices
    if size(viRow,2)==1, viRow=viRow'; end %row
    viSpk = int32(spkLim(1):spkLim(end))';
    miRange = bsxfun(@plus, viSpk, int32(viRow));
    miRange = min(max(miRange, 1), numel(vr));
    if isempty(viCol)
        mr = vr(miRange); %2x faster
    else
        mr = vr(miRange, viCol); %matrix passed to save time
    end
end
