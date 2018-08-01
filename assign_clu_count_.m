%--------------------------------------------------------------------------
function S_clu = assign_clu_count_(S_clu, P)
    nRepeat_max = 1000;
    if isempty(P.minClusterSize), P.minClusterSize = 0; end
    if ~isfield(S_clu, 'viClu'), S_clu.viClu = []; end
    if isempty(S_clu.viClu)
        nClu_pre = [];
    else
        nClu_pre = S_clu.nClusters;
    end
    nClu_rm = 0;
    fprintf('assigning clusters, nClu:%d\n', numel(S_clu.icl)); t1=tic;

    if get_set_(P, 'f_assign_site_clu', 0)
        S_clu.viClu = assignCluster_site_(S_clu, get0_());
    end

    % fReassign = 0;
    % min_rho = -inf;
    for iRepeat=1:nRepeat_max % repeat 1000 times max
        %     S_clu.icl(S_clu.rho(S_clu.icl) < min_rho) = [];
        [S_clu.viClu, S_clu.icl] = assignCluster_(S_clu.viClu, S_clu.ordrho, S_clu.nneigh, S_clu.icl);
        %     S_clu.viClu(S_clu.rho < min_rho) = 0; %noise assignment
        %     if isempty(P.minClusterSize), P.minClusterSize = 0; end
        P.minClusterSize = max(get_set_(P, 'minClusterSize', 0), S_clu.trFet_dim(1)*2);
        % http://scikit-learn.org/stable/modules/lda_qda.html

        S_clu = S_clu_refresh_(S_clu);

        % remove clusters unused
        viCluKill = find(S_clu.vnSpk_clu <= P.minClusterSize);
        if isempty(viCluKill), break; end
        S_clu.icl(viCluKill) = [];
        S_clu.viClu=[];
        nClu_rm = nClu_rm + numel(viCluKill);
        if iRepeat==nRepeat_max
            fprintf(2, 'assign_clu_count_: exceeded nRepeat_max=%d\n', nRepeat_max);
        end
    end

    fprintf('\n\ttook %0.1fs. Removed %d clusters having <%d spikes: %d->%d\n', ...
    toc(t1), nClu_rm, P.minClusterSize, nClu_pre, S_clu.nClusters);
end
