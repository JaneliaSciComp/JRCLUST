function hFigs = doSpawnFigures(hCfg)
    %DOSPAWNFIGURES Create the standard cadre of figures
    hFigs = containers.Map();

    fSkipRD = ~ismember('FigRD', hCfg.figList);
    if fSkipRD
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

    hFigs('FigPos')    = doCreateFigure('FigPos', [0 0 .15 .5], ['Unit position; ', hCfg.sessionName], 1, 0);
    hFigs('FigMap')    = doCreateFigure('FigMap', [0 .5 .15 .5], ['Probe map; ', hCfg.sessionName], 1, 0);
    hFigs('FigWav')    = doCreateFigure('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', hCfg.sessionName], 0, 1);
    hFigs('FigTime')   = doCreateFigure('FigTime', ftShape, ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', hCfg.sessionName]);
    hFigs('FigProj')   = doCreateFigure('FigProj', [.5 .2 .35 .5], ['Feature projection: ', hCfg.sessionName]);
    hFigs('FigSim')    = doCreateFigure('FigSim', [.5 .7 .35 .3], ['Waveform-based similarity score (click): ', hCfg.sessionName]);
    hFigs('FigHist')   = doCreateFigure('FigHist', fhShape, ['ISI Histogram: ', hCfg.sessionName]);
    hFigs('FigISI')    = doCreateFigure('FigISI', fiShape, ['Return map: ', hCfg.sessionName]);
    hFigs('FigCorr')   = doCreateFigure('FigCorr', fcShape, ['Time correlation: ', hCfg.sessionName]);
    if ~fSkipRD
        hFigs('FigRD') = doCreateFigure('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', hCfg.sessionName]);
    end
end
