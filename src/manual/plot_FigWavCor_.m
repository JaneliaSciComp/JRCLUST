%--------------------------------------------------------------------------
function [hFig, S_fig] = plot_FigWavCor_(S0, figView)
    if nargin < 1
        S0 = get(0, 'UserData');
    end
    if nargin < 2
        figView = 'wavecor';
    end

    S_clu = S0.S_clu;
    P = S0.P;
    [hFig, S_fig] = get_fig_cache_('FigWavCor');

    figure_wait_(1, hFig);
    nClu = S_clu.nClu;

    % Plot
    if isempty(S_fig)
        S_fig.hAx = axes_new_(hFig);
        set(S_fig.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
        set(S_fig.hAx, {'XTick', 'YTick'}, {1:nClu, 1:nClu});

        if get_set_(P, 'fImportKsort', 0)
            S_fig.figView = 'simscore';
        else
            S_fig.figView = 'wavecor';
        end

        axis_(S_fig.hAx, [0 nClu 0 nClu] + .5);
        axis(S_fig.hAx, 'xy');
        grid(S_fig.hAx, 'on');
        xlabel(S_fig.hAx, 'Clu#');
        ylabel(S_fig.hAx, 'Clu#');
        if strcmp(S_fig.figView, 'simscore')
            S_fig.hImWavCor = imagesc(S_clu.mrSim_clu, P.corrLim); %clears title and current figure
        else
            S_fig.hImWavCor = imagesc(S_clu.mrWavCor, P.corrLim); %clears title and current figure
        end
        S_fig.hCursorV = line([1 1], [.5 nClu + .5], 'Color', [0 0 0], 'LineWidth', 1.5);
        S_fig.hCursorH = line([.5 nClu + .5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);
        colorbar(S_fig.hAx);

        if get_set_(P, 'fImportKsort', 0)
            S_fig.vcTitle = '[S]plit; [M]erge; [D]elete; [K]iloSort sim score; [W]aveform corr';
        else
            S_fig.vcTitle = '[S]plit; [M]erge; [D]elete';
        end

        set(hFig, 'KeyPressFcn', @keyPressFcn_FigWavCor_);
        mouse_figure(hFig, S_fig.hAx, @button_FigWavCor_);
        S_fig.hDiag = plotDiag_([0, nClu, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        S_fig.figView = figView;
        if strcmp(S_fig.figView, 'simscore')
            set(S_fig.hImWavCor, 'CData', S_clu.mrSim_clu);
            set(hFig, 'Name', ['KiloSort cluster similarity score (click): ', P.vcFile_prm], 'NumberTitle', 'off', 'Color', 'w');
        else
            set(S_fig.hImWavCor, 'CData', S_clu.mrWavCor);
            set(hFig, 'Name', ['Waveform correlation (click): ', P.vcFile_prm], 'NumberTitle', 'off', 'Color', 'w');
        end

        set(S_fig.hAx, {'XTick', 'YTick'}, {1:nClu, 1:nClu});
        [diagX, diagY] = plotDiag__([0, nClu, .5]);
        set(S_fig.hDiag, 'XData', diagX, 'YData', diagY);
    end

    % output
    set(hFig, 'UserData', S_fig);
    figure_wait_(0, hFig);
end %func
