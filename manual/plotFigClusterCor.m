%--------------------------------------------------------------------------
function [hFig, figData] = plotFigClusterCor(S0, figView)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    if nargin < 2
        figView = 'wavecor';
    end

    fImportKilosort = getOr(S0.P, 'fImportKilosort', 0);

    S_clu = S0.S_clu;
    P = S0.P;
    [hFig, figData] = getCachedFig('FigClusterCor');

    figure_wait_(1, hFig);
    nClusters = S_clu.nClusters;

    % Plot
    if isempty(figData)
        figData.hAx = axes_new_(hFig);
        set(figData.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
        set(figData.hAx, {'XTick', 'YTick'}, {1:nClusters, 1:nClusters});

        if fImportKilosort
            figData.figView = 'simscore';
        else
            figData.figView = 'wavecor';
        end

        axis_(figData.hAx, [0 nClusters 0 nClusters] + .5);
        axis(figData.hAx, 'xy');
        grid(figData.hAx, 'on');
        xlabel(figData.hAx, 'Clu#');
        ylabel(figData.hAx, 'Clu#');
        if strcmp(figData.figView, 'simscore')
            figData.hImClusterCor = imagesc(S_clu.simScore, P.corrLim); %clears title and current figure
        else
            figData.hImClusterCor = imagesc(S_clu.mrWavCor, P.corrLim); %clears title and current figure
        end
        figData.hCursorV = line([1 1], [.5 nClusters+.5], 'Color', [0 0 0], 'LineWidth', 1.5);
        figData.hCursorH = line([.5 nClusters+.5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);
        colorbar(figData.hAx);

        if fImportKilosort
            figData.vcTitle = '[S]plit; [M]erge; [D]elete; [T]emplate similarity; [W]aveform correlation';
        else
            figData.vcTitle = '[S]plit; [M]erge; [D]elete';
        end

        set(hFig, 'KeyPressFcn', @keyPressFigClusterCor);
        mouse_figure(hFig, figData.hAx, @buttonFigClusterCor);
        figData.hDiag = plotDiag_([0, nClusters, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        fSwapView = ~strcmp(figData.figView, figView);

        figData.figView = figView;
        if strcmp(figData.figView, 'simscore')
            set(figData.hImClusterCor, 'CData', S_clu.simScore);
            set(hFig, 'Name', ['KiloSort cluster similarity score (click): ', P.paramFile], 'NumberTitle', 'off', 'Color', 'w');
        else
            set(figData.hImClusterCor, 'CData', S_clu.mrWavCor);
            set(hFig, 'Name', ['Waveform correlation (click): ', P.paramFile], 'NumberTitle', 'off', 'Color', 'w');
        end

        if ~fSwapView
            set(figData.hCursorV, 'xdata', [1 1], 'ydata', [.5 nClusters+.5]);
            set(figData.hCursorH, 'xdata', .5+[0 nClusters], 'ydata', [1 1]);
        end
    end

    % output
    set(hFig, 'UserData', figData);
    figure_wait_(0, hFig);
end %func
