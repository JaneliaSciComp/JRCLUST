%--------------------------------------------------------------------------
function FigClusterCor_update_(S0)
    if nargin<1, S0 = get(0, 'UserData'); end

    [hFig, S_fig] = getCachedFig('FigClusterCor');
    set(S_fig.hImWavCor, 'CData', S0.S_clu.mrWavCor);
    % plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5, 'Parent', S_fig.); %plot in one scoop
    nClu = size(S0.S_clu.mrWavCor, 1);
    [vrX, vrY] = plotDiag__([0, nClu, .5]);
    set(S_fig.hDiag, 'XData', vrX, 'YData', vrY);
end %func
