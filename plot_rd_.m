%--------------------------------------------------------------------------
function hFig = plot_rd_(P, S0)
    fLog_delta = 0;
    max_delta = 5;
    if nargin<2, S0 = []; end
    if isempty(S0), S0 = load0_(P); end
    S_clu = S0.S_clu;
    % S_clu = get0_('S_clu');
    if fLog_delta
        vrY_plot = log10(S_clu.delta);
        vcY_label = 'log10 Delta';
        vrY_plot(vrY_plot>log10(max_delta)) = nan;
    else
        vrY_plot = (S_clu.delta);
        vcY_label = 'Delta';
        vrY_plot(vrY_plot>max_delta) = nan;
    end
    vrX_plot = S_clu.rho;
    vrX_plot1 = log10(S_clu.rho);
    [~, ~, vrY_plot1] = detrend_local_(S_clu, P, 0);
    vrY_plot1 = log10(vrY_plot1);

    % Plot
    vhAx = [];
    hFig = createFigure('Fig_RD1', [0 0 .5 1], P.paramFile, 1, 1);
    vhAx(1) = subplot(211); hold on;
    plot(vrX_plot1, vrY_plot, '.', 'Color', repmat(.5,1,3));
    plot(vrX_plot1(S_clu.icl), vrY_plot(S_clu.icl), 'ro');
    plot(P.log10RhoCutoff*[1,1], [0, max_delta], 'r-');
    xylabel_(gca, 'log10 Rho', vcY_label, sprintf('lin-log plot, nClu: %d', S_clu.nClusters));
    axis_([-4, -.5, 0, max_delta]); grid on;

    vhAx(2) = subplot(212); hold on;
    plot(vrX_plot1, vrY_plot1, '.', 'Color', repmat(.5,1,3));
    plot(vrX_plot1(S_clu.icl), vrY_plot1(S_clu.icl), 'ro');
    plot(P.log10RhoCutoff*[1,1], [-.5, 2], 'r-', [-4, -.5], P.log10DeltaCutoff*[1,1], 'r-');
    xylabel_(gca, 'log10 Rho', 'log10 Delta (detrended)', ...
    sprintf('detrended log-log plot (P.log10RhoCutoff=%0.3f; P.log10DeltaCutoff=%0.3f; nClu:%d)', ...
    P.log10RhoCutoff, P.log10DeltaCutoff, S_clu.nClusters));
    axis_([-4, -.5, -.5, 2]); grid on;

    linkaxes(vhAx, 'x');
end %func
