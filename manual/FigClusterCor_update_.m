%--------------------------------------------------------------------------
function FigClusterCor_update_(S0)
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    [hFig, S_fig] = getCachedFig('FigClusterCor');
    set(S_fig.hImWavCor, 'CData', S0.S_clu.mrWavCor);

    nClu = size(S0.S_clu.mrWavCor, 1);
    [vrX, vrY] = getGridDiagonal([0, nClu, .5]);
    set(S_fig.hDiag, 'XData', vrX, 'YData', vrY);
end %func
