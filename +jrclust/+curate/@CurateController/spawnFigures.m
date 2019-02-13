function spawnFigures(obj)
    %SPAWNFIGURES Create the standard cadre of figures
    hFigs_ = containers.Map();

    skipRD = ~ismember('FigRD', obj.hCfg.figList); % skip the rho-delta plot
    if skipRD % expand figTime et al. to take its place
        ftShape = [.15 0 .85 .2];
        fcShape = [.85 .2 .15 .27];
        fiShape = [.85 .47 .15 .26];
        fhShape = [.85 .73 .15 .27];
    else
        ftShape = [.15 0 .7 .2];
        fcShape = [.85 .25 .15 .25];
        fiShape = [.85 .5 .15 .25];
        fhShape = [.85 .75 .15 .25];
    end

    hFigs_('FigPos')    = createFig('FigPos', [0 0 .15 .5], ['Unit position; ', obj.hCfg.sessionName], 1, 0);
    hFigs_('FigMap')    = createFig('FigMap', [0 .5 .15 .5], ['Probe map; ', obj.hCfg.sessionName], 1, 0);
    hFigs_('FigWav')    = createFig('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', obj.hCfg.sessionName], 0, 1);
    hFigs_('FigTime')   = createFig('FigTime', ftShape, ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', obj.hCfg.sessionName]);
    hFigs_('FigProj')   = createFig('FigProj', [.5 .2 .35 .5], ['Feature projection: ', obj.hCfg.sessionName]);
    hFigs_('FigSim')    = createFig('FigSim', [.5 .7 .35 .3], ['Waveform-based similarity score (click): ', obj.hCfg.sessionName]);
    hFigs_('FigHist')   = createFig('FigHist', fhShape, ['ISI Histogram: ', obj.hCfg.sessionName]);
    hFigs_('FigISI')    = createFig('FigISI', fiShape, ['Return map: ', obj.hCfg.sessionName]);
    hFigs_('FigCorr')   = createFig('FigCorr', fcShape, ['Time correlation: ', obj.hCfg.sessionName]);
    if ~skipRD
        hFigs_('FigRD') = createFig('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', obj.hCfg.sessionName]);
    end

    obj.hFigs = hFigs_;
end

%% LOCALFUNCTIONS
function hFig = createFig(figTag, figPos, figName, figToolbar, figMenubar)
    %CREATEFIG Create a Figure object
    if nargin < 2
        figPos = [];
    end
    if nargin < 3
        figName = '';
    end
    if nargin < 4
        figToolbar = 0;
    end
    if nargin < 5
        figMenubar = 0;
    end

    hFig = jrclust.views.Figure(figTag, figPos, figName, figToolbar, figMenubar);
end
