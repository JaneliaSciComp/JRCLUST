%--------------------------------------------------------------------------
function S_clu = plotFigRD(S_clu, P)
    hFig= getCachedFig('FigRD');
    figure(hFig);
    clf;

    % TODO: cS_clu_shank not referred to anywhere else; remove
    if isfield(S_clu, 'cS_clu_shank')
        cellfun(@(S_clu1) plotFigRD(S_clu1, P), S_clu.cS_clu_shank);
        return;
    end

    if isempty(P.delta1_cut)
        P.delta1_cut = S_clu.P.delta1_cut;
    end

    if isempty(P.rho_cut)
        P.rho_cut = S_clu.P.rho_cut;
    end

    if isempty(P.min_count)
        P.min_count = S_clu.P.min_count;
    end

    if ~isfield(P, 'vcDetrend_postclu')
        P.vcDetrend_postclu = 'none';
    end

    switch P.vcDetrend_postclu
        case 'none'
            icl = find(S_clu.rho(:) > 10^(P.rho_cut) & S_clu.delta(:) > 10^(P.delta1_cut));
            x = log10_(S_clu.rho(:));
            y = log10_(S_clu.delta(:));
            fDetrend = 0;
        case 'global'
            [icl, x, y] = detrend_local_(S_clu, P, 0);
            vl_nan = y<=0;
            y = log10_(y);
            fDetrend = 1;
        case 'local'
            [icl, x, y] = detrend_local_(S_clu, P, 1);
            y = log10_(y);
            fDetrend = 1;
    end % switch

    hold on; plot(x, y, '.');
    axis tight;
    axis_([-4 -.5 -1 2])
    set(gcf,'color','w');
    set(gcf, 'UserData', struct('x', x, 'y', y)); grid on;
    set(gca,'XScale','linear', 'YScale', 'linear');
    plot(P.rho_cut*[1 1], get(gca,'YLim'), 'r--', get(gca,'XLim'), P.delta1_cut*[1, 1], 'r--');
    xlabel('log10 rho'); ylabel(sprintf('log10 delta (detrend=%d)', fDetrend));

    % label clusters
    if isfield(S_clu, 'icl')
        icl = S_clu.icl; % do not overwrite
    end
    x_icl = double(x(icl));
    y_icl = double(y(icl));
    % if P.fLabelClu
    %     arrayfun(@(i)text(x_icl(i), y_icl(i), sprintf('%dn%d',i,S_clu.vnSpk_clu(i)), 'VerticalAlignment', 'bottom'), 1:numel(icl));
    % end
    hold on;
    plot(x_icl, y_icl, 'r.');
    grid on;
    % nClu = numel(unique(S_clu.viClu(S_clu.viClu>0))); %numel(icl)
    title_(sprintf('rho-cut:%f, delta-cut:%f', P.rho_cut, P.delta1_cut));
    drawnow;
end %func
