%--------------------------------------------------------------------------
function S0 = figures_manual_(P)
    % 'iFig', [], 'name', '', 'pos', [], 'fToolbar', 0, 'fMenubar', 0);

    csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigCorr', 'FigIsi', 'FigHist'};

    if get_set_(P, 'fImportKsort', 0)
        ftShape = [.15 0 .85 .2];
        fcShape = [.85 .2 .15 .27];
        fiShape = [.85 .47 .15 .26];
        fhShape = [.85 .73 .15 .27];
        fwcTitle = ['KiloSort cluster similarity score (click): ', P.vcFile_prm];
    else
        ftShape = [.15 0 .7 .2];
        fcShape = [.85 .25 .15 .25];
        fiShape = [.85 .5 .15 .25];
        fhShape = [.85 .75 .15 .25];
        fwcTitle = ['Waveform correlation (click): ', P.vcFile_prm];

        % rho-delta plot
        csFig{end+1} = 'FigRD';
        create_figure_('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.vcFile_prm]);
    end

    create_figure_('FigPos', [0 0 .15 .5], ['Unit position; ', P.vcFile_prm], 1, 0);

    create_figure_('FigMap', [0 .5 .15 .5], ['Probe map; ', P.vcFile_prm], 1, 0);

    create_figure_('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', P.vcFile_prm], 0, 1);

    create_figure_('FigTime', ftShape, ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

    create_figure_('FigProj', [.5 .2 .35 .5], ['Feature projection: ', P.vcFile_prm]);

    create_figure_('FigWavCor', [.5 .7 .35 .3], fwcTitle);

    create_figure_('FigHist', fhShape, ['ISI Histogram: ', P.vcFile_prm]);

    create_figure_('FigIsi', fiShape, ['Return map: ', P.vcFile_prm]);

    create_figure_('FigCorr', fcShape, ['Time correlation: ', P.vcFile_prm]);

    cvrFigPos0 = cellfun(@(vc) get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);

    S0 = set0_(cvrFigPos0, csFig);
end %func
