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
    axOffset = 0.8;
    axLen = 1/nStims;

    if ~jrclust.utils.isvalid(hFigTrial1)
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

    if ~jrclust.utils.isvalid(hFigTrial2)
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

    for iStim = 1:n
        cla(vhAx1(iStim));
        cla(vhAx2(iStim));
        vrTime_trial = trialTimes{iStim}; %(:,1);
        nTrials = numel(vrTime_trial);
        viTime_clu1 = S_clu_time_(hClust, iCluster);
        plot_raster_clu_(viTime_clu1, vrTime_trial, hCfg, vhAx1(iStim));
        plot_psth_clu_(viTime_clu1, vrTime_trial, hCfg, vhAx2(iStim), hFigTrial.figData.color);
        title(vhAx2(iStim), sprintf('Cluster %d; %d trials', iCluster, nTrials));
    end
    %     offset = offset + nTrials;
    if numel(vhAx1)>2
        set(vhAx1(2:end),'xticklabel',{});
        for ax = vhAx1(2:end)
            xlabel(ax, '')
        end
    end % end
end %func

function plot_raster_clu_(viTime_clu, vrTime_trial, P, hAx)
    if nargin<4, hAx=gca; end

    trialLength = diff(P.tlim_psth); % seconds
    nTrials = numel(vrTime_trial);
    spikeTimes = cell(nTrials, 1);
    t0 = -P.tlim_psth(1);
    for iTrial = 1:nTrials
        rTime_trial1 = vrTime_trial(iTrial);
        vrTime_lim1 = rTime_trial1 + P.tlim_psth;
        vrTime_clu1 = double(viTime_clu) / P.sRateHz;
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
end

function plot_psth_clu_(viTime_clu, vrTime_trial, P, hAx, vcColor)
    if nargin<4, hAx=gca; end
    if nargin<5, vcColor = 'k'; end

    tbin = P.tbin_psth;
    nbin = round(tbin * P.sRateHz);
    nlim = round(P.tlim_psth/tbin);
    viTime_Trial = round(vrTime_trial / tbin);

    vlTime1=zeros(0);
    vlTime1(ceil(double(viTime_clu)/nbin))=1;
    mr1 = vr2mr2_(double(vlTime1), viTime_Trial, nlim);
    vnRate = mean(mr1,2) / tbin;
    vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;
    bar(hAx, vrTimePlot, vnRate, 1, 'EdgeColor', 'none', 'FaceColor', vcColor);
    vrXTick = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
    set(hAx, 'XTick', vrXTick, 'XTickLabel', []);
    grid(hAx, 'on');
    hold(hAx, 'on');
    plot(hAx, [0 0], get(hAx,'YLim'), 'r-');
    ylabel(hAx, 'Rate (Hz)');
    xlim_(hAx, P.tlim_psth);
end
