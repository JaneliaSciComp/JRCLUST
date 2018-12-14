%--------------------------------------------------------------------------
function S_clu = postCluster_(S_clu, P)
    % fSortClu = 0; %debug
    if isfield(S_clu, 'viClu')
        S_clu = rmfield(S_clu, 'viClu');
    end
    if isempty(S_clu), return; end

    switch lower(P.vcDetrend_postclu)
        case {'hidehiko', 'hide'}
        S_clu.icl = select_icl_with_slope_(S_clu, P.delta1_cut);
        case 'local' %
        S_clu.icl = detrend_local_(S_clu, P, 1);
        case 'none'
        S_clu.icl = find(S_clu.rho(:) > 10^(P.rho_cut) & S_clu.delta(:) > 10^(P.delta1_cut));
        case 'global'
        S_clu.icl = detrend_local_(S_clu, P, 0);
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
