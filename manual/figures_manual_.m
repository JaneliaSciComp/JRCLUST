%--------------------------------------------------------------------------
function S0 = figures_manual_(P)
    % 'iFig', [], 'name', '', 'pos', [], 'fToolbar', 0, 'fMenubar', 0);
    create_figure_('FigPos', [0 0 .15 .5], ['Unit position; ', P.vcFile_prm], 1, 0);
    create_figure_('FigMap', [0 .5 .15 .5], ['Probe map; ', P.vcFile_prm], 1, 0);

    create_figure_('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', P.vcFile_prm], 0, 1);
    create_figure_('FigTime', [.15 0 .7 .2], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

    create_figure_('FigProj', [.5 .2 .35 .5], ['Feature projection: ', P.vcFile_prm]);
    create_figure_('FigWavCor', [.5 .7 .35 .3], ['Waveform correlation (click): ', P.vcFile_prm]);

    create_figure_('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.vcFile_prm]);
    create_figure_('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.vcFile_prm]);
    create_figure_('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.vcFile_prm]);
    create_figure_('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.vcFile_prm]);

    csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
    cvrFigPos0 = cellfun(@(vc)get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);
    S0 = set0_(cvrFigPos0, csFig);
end %func
