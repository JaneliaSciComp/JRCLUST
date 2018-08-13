%--------------------------------------------------------------------------
function S_clu = plotFigRD(S_clu, P)
    hFig= getCachedFig('FigRD');
    set(hFig, 'DefaultTextInterpreter', 'LaTeX');
    figure(hFig);
    clf;

    % TODO: cS_clu_shank not referred to anywhere else; remove?
    % if isfield(S_clu, 'cS_clu_shank')
    %     cellfun(@(S_clu1) plotFigRD(S_clu1, P), S_clu.cS_clu_shank);
    %     return;
    % end

    if ~isfield(P, 'log10DeltaCutoff') || isempty(P.log10DeltaCutoff)
        P.log10DeltaCutoff = S_clu.P.log10DeltaCutoff;
    end

    if ~isfield(P, 'log10RhoCutoff') || isempty(P.log10RhoCutoff)
        P.log10RhoCutoff = S_clu.P.log10RhoCutoff;
    end

    if ~isfield(P, 'minClusterSize') || isempty(P.minClusterSize)
        P.minClusterSize = S_clu.P.minClusterSize;
    end

    if ~isfield(P, 'vcDetrend_postclu')
        P.vcDetrend_postclu = 'none';
    end

    clusterCenters = S_clu.icl;

    switch P.vcDetrend_postclu
        case 'none'
            x = nanLog10(S_clu.rho(:));
            y = nanLog10(S_clu.delta(:));
            fDetrend = 0;
        case 'global'
            [clusterCenters, x, y] = detrend_local_(S_clu, P, 0);
            y = nanLog10(y);
            fDetrend = 1;
        case 'local'
            [clusterCenters, x, y] = detrend_local_(S_clu, P, 1);
            y = nanLog10(y);
            fDetrend = 1;
    end % switch

    hold on;
    plot(x, y, '.');
    axis tight;

    axis_([-4 -.5 -1 2])
    set(gcf,'color','w');
    set(gcf, 'UserData', struct('x', x, 'y', y)); grid on;
    set(gca, 'XScale', 'linear', 'YScale', 'linear');

    plot(P.log10RhoCutoff * [1 1], get(gca, 'YLim'), 'r--', ...
         get(gca, 'XLim'), P.log10DeltaCutoff * [1, 1], 'r--');
    xlabel('$\log_{10} \rho$');
    ylab = '$\log_{10} \delta$';
    if fDetrend
        ylab = [ylab '( detrended)'];
    end
    ylabel(ylab);

    % plot cluster centers
    ccX = double(x(clusterCenters));
    ccY = double(y(clusterCenters));

    hold on;
    plot(ccX, ccY, 'r.');
    grid on;

    title_(sprintf('rho-cut: %f, delta-cut: %f', P.log10RhoCutoff, P.log10DeltaCutoff));

    drawnow;
end %func
