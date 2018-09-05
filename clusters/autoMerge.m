%--------------------------------------------------------------------------
function autoMerge(S0)
    % UI function; merge clusters with correlation above a user-input threshold

    if nargin < 1 || isempty(S0)
        S0 = get(0, 'UserData');
    end
    S_clu = S0.S_clu;
    P = S0.P;

    % get user input
    response = inputdlg_('Waveform correlation threshold (0-1):', 'Auto-merge based on waveform threshold', 1, {num2str(P.maxWavCor)});
    if isempty(response)
        return;
    end
    maxWavCor = str2double(response{1});

    % sanity check input value
    if isnan(maxWavCor) || maxWavCor < 0 || maxWavCor > 1
        msgbox_('Invalid criterion.');
        return;
    end

    figure_wait_(1);
    drawnow;

    nClustersBefore = S_clu.nClusters;
    P.maxWavCor = maxWavCor;

    S_clu = post_merge_wav_(S_clu, P.nRepeat_merge, P);
    % [S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, [], S0));
    setUserData(S_clu);

    S0 = gui_update_();
    figure_wait_(0);

    dialogAssert(clusterDataConsistent(S_clu), 'Cluster number is inconsistent after deleting');
    nClu_merge = nClustersBefore - S_clu.nClusters;
    msgbox_(sprintf('Merged %d clusters >%0.2f maxWavCor.', nClu_merge, maxWavCor));
    save_log_(sprintf('merge-auto <%0.2f maxWavCor', maxWavCor), S0);
end % function
