%--------------------------------------------------------------------------
% 10/27/17 JJJ: Delete multiple clusters
function S_clu = delete_clu_(S_clu, viClu_delete)
    % sets the cluster to zero
    nClu_prev = S_clu.nClu;
    viClu_keep = setdiff(1:nClu_prev, viClu_delete);
    try
        S_clu = S_clu_select_(S_clu, viClu_keep); % remap all
    catch
        disp('err');
    end

    iClu_del = min(S_clu.viClu) - 1;
    if iClu_del==0, iClu_del = -1; end
    vlDelete_spk = ismember(S_clu.viClu, viClu_delete);
    S_clu.viClu(vlDelete_spk) = iClu_del;
    nClu_new = numel(viClu_keep);

    vlMap = S_clu.viClu > 0;
    viMap = zeros(1, nClu_prev);
    viMap(viClu_keep) = 1:nClu_new;
    S_clu.viClu(vlMap) = viMap(S_clu.viClu(vlMap));
    S_clu.nClu = nClu_new;
    % update viClu
    % if viClu_delete < max(S_clu.viClu)
    %     viUpdate = find(S_clu.viClu>viClu_delete);
    %     S_clu.viClu(viUpdate) = S_clu.viClu(viUpdate) - 1;
    % end
    % for iClu3 = viClu_delete+1:S_clu.nClu % update cluster chain info
    %     S_clu = S_clu_update_note_(S_clu, iClu3, get_next_clu_(S_clu, iClu3) - 1);
    % end
    assert_(S_clu_valid_(S_clu), 'Cluster number is inconsistent after deleting');
end %func
