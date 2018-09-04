%--------------------------------------------------------------------------
function S_clu = postCluster_(S_clu, P)
    % fSortClu = 0; %debug
    if isfield(S_clu, 'spikeClusters')
        S_clu = rmfield(S_clu, 'spikeClusters');
    end
    if isempty(S_clu), return; end

    switch lower(P.vcDetrend_postclu)
        %         case {'hidehiko', 'hide'}
        %             S_clu.clusterCenters = selec_rho_delta_with_slope(S_clu, P.log10DeltaCutoff);
        case 'local' %
            S_clu.clusterCenters = detrend_local_(S_clu, P, 1);
        case 'none'
            S_clu.clusterCenters = find(S_clu.rho(:) > 10^(P.log10RhoCutoff) & S_clu.delta(:) > 10^(P.log10DeltaCutoff));
        case 'global'
            S_clu.clusterCenters = detrend_local_(S_clu, P, 0);
        case 'logz'
            S_clu.clusterCenters = log_ztran_(S_clu.rho, S_clu.delta, P.log10RhoCutoff, 4+P.log10DeltaCutoff);
        otherwise
            fprintf(2, 'postCluster_: vcDetrend_postclu = ''%s''; not supported.\n', P.vcDetrend_postclu);
    end
    % end

    % Update P
    S_clu.P.minClusterSize = P.minClusterSize;
    S_clu.P.log10DeltaCutoff = P.log10DeltaCutoff;
    S_clu.P.log10RhoCutoff = P.log10RhoCutoff;
    S_clu.spikeClusters = [];
    S_clu = assign_clu_count_(S_clu, P); % enforce min count algorithm

    % debug output
    if nargout==0
        vrXp = log10(S_clu.rho);
        vrYp = log10(S_clu.delta);
        figure; hold on;
        plot(vrXp, vrYp, 'b.', vrXp(S_clu.clusterCenters), vrYp(S_clu.clusterCenters), 'ro');
        plot(get(gca, 'XLim'), P.log10DeltaCutoff*[1 1], 'r-');
        plot(P.log10RhoCutoff*[1 1], get(gca, 'YLim'), 'r-');
        xlabel('log10 rho');
        ylabel('log10 delta');
    end
end % function
