%--------------------------------------------------------------------------
function S_clu = merge_clu_(S_clu, iClu1, iClu2, P)
    if iClu1>iClu2, [iClu1, iClu2] = swap_(iClu1, iClu2); end

    S_clu = merge_clu_pair_(S_clu, iClu1, iClu2);
    S_clu = S_clu_refrac_(S_clu, P, iClu1); % remove refrac
    S_clu = S_clu_update_(S_clu, iClu1, P);
    S_clu = delete_clu_(S_clu, iClu2);
    % S_clu = S_clu_remove_empty_(S_clu);
    dialogAssert(clusterDataConsistent(S_clu), 'Cluster number is inconsistent after merging');
    fprintf('%s [W] merging Clu %d and %d\n', datestr(now, 'HH:MM:SS'), iClu1, iClu2);
end % function
