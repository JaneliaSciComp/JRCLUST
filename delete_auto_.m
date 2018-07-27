%--------------------------------------------------------------------------
function delete_auto_()
    % SNR based delete functionality
    % Ask SNR
    S0 = get(0, 'UserData');
    [S_clu, P] = get0_('S_clu', 'P');
    hFig = create_figure_('', [.5 .7 .35 .3], ['Delete Auto: ', P.vcFile]);

    % Ask user which clusters to delete
    plot(S_clu.vrSnr_clu(:), S_clu.vnSpk_clu(:), '.'); % show cluster SNR and spike count
    xlabel('Unit SNR'); ylabel('# spikes/unit'); grid on;
    set(gca,'YScale','log');
    % snr_thresh = inputdlg_num_('SNR threshold: ', 'Auto-deletion based on SNR', 10); % also ask about # spikes/unit (or firing rate) @TODO
    csAns = inputdlg_({'Min Unit SNR:', 'Max Unit SNR:', 'Minimum # spikes/unit'}, 'Auto-deletion based on SNR', 1, {'5', 'inf', '0'}); % also ask about # spikes/unit (or firing rate) @TODO
    close(hFig);

    % parse user input
    if isempty(csAns), return; end
    snr_min_thresh = str2double(csAns{1});
    snr_max_thresh = str2double(csAns{2});
    count_thresh = round(str2double(csAns{3}));
    if any(isnan([snr_min_thresh, snr_max_thresh, count_thresh]))
        msgbox_('Invalid criteria.'); return;
    end
    viClu_delete = find(S_clu.vrSnr_clu(:) < snr_min_thresh | S_clu.vnSpk_clu(:) < count_thresh | S_clu.vrSnr_clu(:) > snr_max_thresh);
    if isempty(viClu_delete), msgbox_('No clusters deleted.'); return; end
    if numel(viClu_delete) >= S_clu.nClu, msgbox_('Cannot delete all clusters.'); return; end

    % Auto delete
    figure_wait_(1); drawnow;
    S_clu = delete_clu_(S_clu, viClu_delete);
    set0_(S_clu);
    S0 = gui_update_();
    figure_wait_(0);

    msgbox_(sprintf('Deleted %d clusters <%0.1f SNR or <%d spikes/unit.', numel(viClu_delete), snr_min_thresh, count_thresh));
    save_log_(sprintf('delete-auto <%0.1f SNR or <%d spikes/unit', snr_min_thresh, count_thresh), S0);
end %func
