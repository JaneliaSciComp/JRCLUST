%--------------------------------------------------------------------------
function S0 = plot_FigRD_(S0)    
    P = S0.P;
    S_clu = S0.S_clu;

    if ~isfield(S0, 'hFigRD')
        return;
    end
    
    S0.hFigRD.clf();

    % chickenshit
%     if isfield(S_clu, 'cS_clu_shank')
%         cellfun(@(S_clu1) plot_FigRD_(S_clu1, P), S_clu.cS_clu_shank);
%         return;
%     end

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
            x = nanlog10(S_clu.rho(:));
            y = nanlog10(S_clu.delta(:));
            fDetrend = 0;

        case 'global'
            [icl, x, y] = detrend_local_(S_clu, P, 0);
            y = nanlog10(y);
            fDetrend = 1;

        case 'local'
            [icl, x, y] = detrend_local_(S_clu, P, 1);
            y = nanlog10(y);
            fDetrend = 1;
    end
    
    S0.hFigRD.hold('on');
    S0.hFigRD.plot(x, y, '.');
    S0.hFigRD.axis('tight');
    S0.hFigRD.axis([-4 -.5 -1 2]);
    S0.hFigRD.figSet('Color', 'w'); % redundant?
    S0.hFigRD.axSet('XScale', 'linear', 'YScale', 'linear');
    S0.hFigRD.plot(P.rho_cut*[1 1], S0.hFigRD.axGet('YLim'), 'r--', ...
                   S0.hFigRD.axGet(gca,'XLim'), P.delta1_cut*[1, 1], 'r--');
    S0.hFigRD.grid('on');

    S0.hFigRD.figMetadata = struct('x', x, 'y', y);
    S0.hFigRD.xlabel('log10 rho');
    S0.hFigRD.ylabel(sprintf('log10 delta (detrend=%d)', fDetrend));

    % label clusters
    if isfield(S_clu, 'icl')
        icl = S_clu.icl; % do not overwrite
    end
    x_icl = double(x(icl));
    y_icl = double(y(icl));
    S0.hFigRD.plot(x_icl, y_icl, 'r.');

    S0.hFigRD.title(sprintf('rho-cut:%f, delta-cut:%f', P.rho_cut, P.delta1_cut));
    drawnow;
    
    % save structs for return
    S0.P = P;
    S0.S_clu = S_clu;
end %func

%% local functions
function z = nanlog10(z)
    z(z <= 0) = nan;
    z = log10(z);
end
