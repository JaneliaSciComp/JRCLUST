%--------------------------------------------------------------------------
% 10/27/17 JJJ: Delete multiple clusters
function S_clu = delete_clu_(S_clu, viClu_delete)
    % sets the cluster to zero
    nClu_prev = S_clu.nClusters;
    viClu_keep = setdiff(1:nClu_prev, viClu_delete);
    try
        S_clu = S_clu_select_(S_clu, viClu_keep); % remap all
    catch
        disp('err');
    end

    iClu_del = min(S_clu.spikeClusters) - 1;
    if iClu_del==0, iClu_del = -1; end
    vlDelete_spk = ismember(S_clu.spikeClusters, viClu_delete);
    S_clu.spikeClusters(vlDelete_spk) = iClu_del;
    nClu_new = numel(viClu_keep);

    vlMap = S_clu.spikeClusters > 0;
    viMap = zeros(1, nClu_prev);
    viMap(viClu_keep) = 1:nClu_new;
    S_clu.spikeClusters(vlMap) = viMap(S_clu.spikeClusters(vlMap));
    S_clu.nClusters = nClu_new;
    % update viClu
    % if viClu_delete < max(S_clu.spikeClusters)
    %     viUpdate = find(S_clu.spikeClusters>viClu_delete);
    %     S_clu.spikeClusters(viUpdate) = S_clu.spikeClusters(viUpdate) - 1;
    % end
    % for iClu3 = viClu_delete+1:S_clu.nClusters % update cluster chain info
    %     S_clu = S_clu_update_note_(S_clu, iClu3, get_next_clu_(S_clu, iClu3) - 1);
    % end
    dialogAssert(S_clu_valid_(S_clu), 'Cluster number is inconsistent after deleting');
end % function
