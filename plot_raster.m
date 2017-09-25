function plot_raster(P, iClu, Sclu)
% import  trial time
% P = loadParam(vcFile_prm);
if nargin < 2, iClu = []; end
if nargin < 3, Sclu = []; end

if isfield(P, 'vcFile_psth'), P.vcFile_trial = P.vcFile_psth; end
%
crTime_trial = loadTrial_(P.vcFile_trial);

%vrTime_trial = loadTrial_(P.vcFile_trial);
if ~iscell(crTime_trial)
    crTime_trial = {crTime_trial};
end
nstims = numel(crTime_trial);
if isempty(crTime_trial), msgbox('Trial file does not exist', 'modal'); return; end

% import cluster time
vcFile_clu = subsFileExt(P.vcFile, '_clu.mat');
if isempty(Sclu), Sclu = loadClu(vcFile_clu); end

% plot lcuster times
Sclu.nClu = max(Sclu.viClu);
% if ~isfield(Sclu, 'cviTime_clu')    
%     Sclu.cviTime_clu = arrayfun(@(iClu)Sclu.viTime(Sclu.viClu==iClu), 1:Sclu.nClu, 'UniformOutput', 0);
% end

hFig = figure(1201); clf;
set(hFig, 'MenuBar', 'none', 'ToolBar', 'none', 'NumberTitle', 'off', 'Color', 'w');
resize_figure(hFig, [.7 .0 .3 .9]);
hTabGroup = uitabgroup(hFig);
% offset = 0;
if isempty(iClu)
    viClu_plot = 1:Sclu.nClu;
else
    viClu_plot = iClu;
end
for iClu = viClu_plot
    htab1 = uitab(hTabGroup, 'Title', sprintf('Clu %d', iClu));    
    axoffset = 0.05;
    axlen = 0.9/nstims;
    hax1 = [];hax2 = [];
    for iStim = 1:nstims
        vrTime_trial = crTime_trial{iStim}(:,1);
        nTrials = numel(vrTime_trial);
        
        % plot raster
        hax1(iStim) = axes('Parent', htab1, 'Position',[.08 axoffset .9 axlen*.68]);
        plot_raster_clu_(Sclu_viTime(Sclu, iClu), vrTime_trial, P);
        
        % plot psth
        hax2(iStim) = axes('Parent', htab1, 'Position',[.08 axoffset + axlen*.68 .9 axlen*.2]);
        plot_psth_clu_(Sclu_viTime(Sclu, iClu), vrTime_trial, P);
        axoffset = axoffset + axlen;
    end
%     offset = offset + nTrials;
    title(sprintf('Cluster %d; %d trials', iClu, nTrials));
    if numel(hax1)>2
        set(hax1(2:end),'xticklabel',{});
        for ax = hax1(2:end)
            axes(ax)
            xlabel('')
        end
        
    end
end
set(gcf, 'Name', P.vcFile);
grid on;
%%
end %func


function plot_psth_clu_(viTime_clu, vrTime_trial, P)
tbin = P.tbin_psth;
nbin = round(tbin * P.sRateHz);
nlim = round(P.tlim_psth/tbin);
viTime_Trial = round(vrTime_trial / tbin);

vlTime1=zeros(0);
vlTime1(ceil(double(viTime_clu)/nbin))=1;
mr1 = vr2mr2(double(vlTime1), viTime_Trial, nlim);
vnRate = mean(mr1,2) / tbin;
vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;

bar(vrTimePlot, vnRate, 1, 'EdgeColor', 'none');

vrXTick = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
set(gca, 'XTick', vrXTick, 'XTickLabel', []);
grid on;
hold on; plot([0 0], get(gca,'YLim'), 'r-');
ylabel('Rate (Hz)');
xlim(P.tlim_psth);
end


function plot_raster_clu_(viTime_clu, vrTime_trial, P)

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
plotSpikeRaster(spikeTimes,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 1], ...
    'LineFormat', struct('LineWidth', 1.5));
ylabel('Trial #')
% title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
vrXTickLabel = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
vrXTick = linspace(0,1,numel(vrXTickLabel));
set(gca, {'XTick', 'XTickLabel'}, {vrXTick, vrXTickLabel});
grid on;
hold on; plot([t0,t0]/trialLength, get(gca,'YLim'), 'r-');
xlabel('Time (s)');
end %func


function vrTime_trial = loadTrial_(vcFile_trial)
% import  trial time
if ~exist(vcFile_trial, 'file'), vrTime_trial = []; return; end

if matchFileExt(vcFile_trial, '.mat')
    Strial = load(vcFile_trial);
    csFields = fieldnames(Strial);
    vrTime_trial = Strial.(csFields{1});
    if isstruct(vrTime_trial)
        vrTime_trial = vrTime_trial.times;
    end
end
end