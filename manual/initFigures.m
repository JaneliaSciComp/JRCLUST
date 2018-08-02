%--------------------------------------------------------------------------
function S0 = initFigures(P) % TODO: create different figures for different algs
    createFigure('FigPos', [0 0 .15 .5], ['Unit position; ', P.paramFile], 1, 0);
    createFigure('FigMap', [0 .5 .15 .5], ['Probe map; ', P.paramFile], 1, 0);

    createFigure('FigTime', [.15 0 .7 .2], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);
    createFigure('FigWav', [.15 .2 .35 .8], ['Averaged waveform: ', P.paramFile], 0, 1);

    createFigure('FigProj', [.5 .2 .35 .5], ['Feature projection: ', P.paramFile]);
    createFigure('FigWavCor', [.5 .7 .35 .3], ['Waveform correlation (click): ', P.paramFile]);

    createFigure('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.paramFile]);
    createFigure('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.paramFile]);
    createFigure('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.paramFile]);
    createFigure('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.paramFile]);

    figTags = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
    figPositions = cellfun(@(c) get(figureByTag(c), 'OuterPosition'), figTags, 'UniformOutput', 0);
    S0 = setUserData(figPositions, figTags);
end % func