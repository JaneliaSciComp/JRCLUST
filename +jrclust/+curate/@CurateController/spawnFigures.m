function spawnFigures(obj)
    %SPAWNFIGURES Create the standard cadre of figures
    hFigs_ = containers.Map();

    skipRD = ~ismember('FigRD', obj.hCfg.figList); % skip the rho-delta plot
    if skipRD % expand figTime et al. to take its place
        ftShape = [0.0015625      0.02963      0.79688      0.19167];
        fcShape = [0.46745      0.69769      0.16406      0.29398];
        fiShape = [.85 .47 .15 .26];
        fhShape = [0.62917      0.69722      0.16823      0.29352];
    else
        ftShape = [0.0013021     0.055556      0.69505      0.18611];
        fcShape = [0.46745      0.69769      0.16406      0.29398];
        fiShape = [.85 .5 .15 .25];
        fhShape = [0.62917      0.69722      0.16823      0.29352];
    end

    hFigs_('FigPos')    = createFig('FigPos', [0.0015625      0.22222      0.12318      0.42778], ['Unit position; ', obj.hCfg.sessionName], 1, 0);
    hFigs_('FigMap')    = createFig('FigMap', [0.0015625      0.65139      0.12344       0.3412], ['Probe map; ', obj.hCfg.sessionName], 1, 0);
    hFigs_('FigWav')    = createFig('FigWav', [0.12266      0.22269      0.34609      0.76944],['Averaged waveform: ', obj.hCfg.sessionName], 0, 1);
    hFigs_('FigTime')   = createFig('FigTime', ftShape, ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', obj.hCfg.sessionName]);
    hFigs_('FigProj')   = createFig('FigProj', [0.46719      0.22407      0.33047      0.47176], ['Feature projection: ', obj.hCfg.sessionName]);
    hFigs_('FigSim')    = createFig('FigSim',[0.79557      0.56759      0.20208      0.42593], ['Waveform-based similarity score (click): ', obj.hCfg.sessionName]);
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
