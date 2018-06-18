%--------------------------------------------------------------------------
function merge_auto_(S0)
    % SNR based delete functionality
    % Ask SNR
    if nargin<1, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [S_clu, P] = deal(S0.S_clu, S0.P);

    % snr_thresh = inputdlg_num_('SNR threshold: ', 'Auto-deletion based on SNR', 10); % also ask about # spikes/unit (or firing rate) @TODO
    csAns = inputdlg_('Waveform correlation threshold (0-1):', 'Auto-merge based on waveform threshold', 1, {num2str(P.maxWavCor)});

    % parse user input
    if isempty(csAns), return; end
    maxWavCor = str2double(csAns{1});
    if isnan(maxWavCor), msgbox_('Invalid criteria.'); return; end

    % Auto delete
    figure_wait_(1); drawnow;
    nClu_prev = S_clu.nClu;
    S_clu = post_merge_wav_(S_clu, P.nRepeat_merge, setfield(P, 'maxWavCor', maxWavCor));
    % [S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
    S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, [], S0));
    set0_(S_clu);
    S0 = gui_update_();
    figure_wait_(0);

    assert_(S_clu_valid_(S_clu), 'Cluster number is inconsistent after deleting');
    nClu_merge = nClu_prev - S_clu.nClu;
    msgbox_(sprintf('Merged %d clusters >%0.2f maxWavCor.', nClu_merge, maxWavCor));
    save_log_(sprintf('merge-auto <%0.2f maxWavCor', maxWavCor), S0);
end %func
