%--------------------------------------------------------------------------
function S0 = initFigures(P) % TODO: create different figures for different algs
    createFigure('FigPos', [0 0 .15 .5], ['Unit position; ', P.paramFile], 1, 0);
    createFigure('FigMap', [0 .5 .15 .5], ['Probe map; ', P.paramFile], 1, 0);

    createFigure('FigTime', [.15 0 .7 .2], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);
    createFigure('FigWav', [.15 .2 .35 .8], ['Averaged waveform: ', P.paramFile], 0, 1);
    createFigure('FigProj', [.5 .2 .35 .5], ['Feature projection: ', P.paramFile]);

    createFigure('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.paramFile]);
    createFigure('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.paramFile]);
    createFigure('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.paramFile]);

    figTags = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigProj', 'FigCorr', 'FigIsi', 'FigHist', 'FigClusterCor'};

    if getOr(P, 'fImportKilosort', 0)
        createFigure('FigClusterCor', [.5 .7 .35 .3], ['KiloSort cluster similarity score (click): ', P.paramFile]);
    else
        createFigure('FigClusterCor', [.5 .7 .35 .3], ['Waveform correlation (click): ', P.paramFile]);
        createFigure('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.paramFile]);
        figTags = {figTags{:}, 'FigRD'};
    end

    figPositions = cellfun(@(c) get(figuresByTag(c), 'OuterPosition'), figTags, 'UniformOutput', 0);

    S0 = setUserData(figPositions, figTags);
end % function
