%--------------------------------------------------------------------------
function [hFig, figData] = plotFigSimScore(S0)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    S_clu = S0.S_clu;
    P = S0.P;
    [hFig, figData] = getCachedFig('FigSimScore');

    figure_wait_(1, hFig);
    nClusters = S_clu.nClusters;

    % Plot
    if isempty(figData)
        figData.hAx = axes_new_(hFig);
        set(figData.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
        set(figData.hAx, {'XTick', 'YTick'}, {1:nClusters, 1:nClusters});
        axis_(figData.hAx, [0 nClusters 0 nClusters] + .5);
        axis(figData.hAx, 'xy');
        grid(figData.hAx, 'on');
        xlabel(figData.hAx, 'Clu#');
        ylabel(figData.hAx, 'Clu#');
        figData.hImSimScore = imagesc(S_clu.simScore, P.corrLim); %clears title and current figure
        figData.hCursorV = line([1 1], [.5 nClusters+.5], 'Color', [0 0 0], 'LineWidth', 1.5);
        figData.hCursorH = line([.5 nClusters+.5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);
        colorbar(figData.hAx);
        figData.vcTitle = '[S]plit; [M]erge; [D]elete';
        set(hFig, 'KeyPressFcn', @keyPressFcn_FigSimScore);
        mouse_figure(hFig, figData.hAx, @button_FigSimScore);
        figData.hDiag = plotDiag_([0, nClusters, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
    else
        set(figData.hImSimScore, 'CData', S_clu.simScore);
        set(figData.hCursorV, 'xdata', [1 1], 'ydata', [.5 nClusters+.5]);
        set(figData.hCursorH, 'xdata', .5+[0 nClusters], 'ydata', [1 1]);
    end
    % output
    set(hFig, 'UserData', figData);
    figure_wait_(0, hFig);
end %func
