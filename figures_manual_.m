%--------------------------------------------------------------------------
function S0 = figures_manual_(P)
    % 'iFig', [], 'name', '', 'pos', [], 'fToolbar', 0, 'fMenubar', 0);
    createFigure('FigPos', [0 0 .15 .5], ['Unit position; ', P.prmFile], 1, 0);
    createFigure('FigMap', [0 .5 .15 .5], ['Probe map; ', P.prmFile], 1, 0);

    createFigure('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', P.prmFile], 0, 1);
    createFigure('FigTime', [.15 0 .7 .2], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

    createFigure('FigProj', [.5 .2 .35 .5], ['Feature projection: ', P.prmFile]);
    createFigure('FigWavCor', [.5 .7 .35 .3], ['Waveform correlation (click): ', P.prmFile]);

    createFigure('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.prmFile]);
    createFigure('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.prmFile]);
    createFigure('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.prmFile]);
    createFigure('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.prmFile]);

    csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
    cvrFigPos0 = cellfun(@(vc)get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);
    S0 = set0_(cvrFigPos0, csFig);
end %func
