%--------------------------------------------------------------------------
function S_clu = postCluster_(S_clu, P)
    % fSortClu = 0; %debug
    if isfield(S_clu, 'viClu')
        S_clu = rmfield(S_clu, 'viClu');
    end
    if isempty(S_clu), return; end

    switch lower(P.vcDetrend_postclu)
        %         case {'hidehiko', 'hide'}
        %             S_clu.icl = selec_rho_delta_with_slope(S_clu, P.delta1_cut);
        case 'local' %
            S0 = get0();
            S_clu.icl = jrclust.clustering.detrendRhoDelta(S_clu, S0.cviSpk_site, true, P);

        case 'none'
            S_clu.icl = find(S_clu.rho(:) > 10^(P.rho_cut) & S_clu.delta(:) > 10^(P.delta1_cut));

        case 'global'
            S0 = get0();
            S_clu.icl = jrclust.clustering.detrendRhoDelta(S_clu, S0.cviSpk_site, false, P);

        case 'logz'
            S_clu.icl = log_ztran_(S_clu.rho, S_clu.delta, P.rho_cut, 4+P.delta1_cut);

        otherwise
            fprintf(2, 'postCluster_: vcDetrend_postclu = ''%s''; not supported.\n', P.vcDetrend_postclu);
    end
    % end

    % Update P
    S_clu.P.min_count = P.min_count;
    S_clu.P.delta1_cut = P.delta1_cut;
    S_clu.P.rho_cut = P.rho_cut;
    S_clu.viClu = [];
    S_clu = assign_clu_count_(S_clu, P); % enforce min count algorithm

    % debug output
    if nargout==0
        vrXp = log10(S_clu.rho);
        vrYp = log10(S_clu.delta);
        figure; hold on;
        plot(vrXp, vrYp, 'b.', vrXp(S_clu.icl), vrYp(S_clu.icl), 'ro');
        plot(get(gca, 'XLim'), P.delta1_cut*[1 1], 'r-');
        plot(P.rho_cut*[1 1], get(gca, 'YLim'), 'r-');
        xlabel('log10 rho');
        ylabel('log10 delta');
    end
end %func
