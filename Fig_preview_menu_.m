%--------------------------------------------------------------------------
% 8/16/17 JJJ: Updated
function Fig_preview_menu_(hFig)
    P = get0_('P');

    set(hFig, 'MenuBar','None');
    mh_file = uimenu(hFig, 'Label','File');
    uimenu(mh_file, 'Label', sprintf('Save to %s', P.paramFile), 'Callback', @(h,e)Fig_preview_save_prm_(hFig));
    uimenu(mh_file, 'Label', '[E]xport to workspace', 'Callback', @(h,e)Fig_preview_export_(hFig));
    uimenu(mh_file, 'Label', 'Save spike detection threshold', 'Callback', @(h,e)Fig_preview_save_threshold_(hFig));

    % Edit menu
    mh_edit = uimenu(hFig, 'Label','Edit');

    uimenu(mh_edit, 'Label', 'Bad site threshold', 'Callback', @(h,e)Fig_preview_site_thresh_(hFig));
    uimenu(mh_edit, 'Label', 'Spike detection threshold', 'Callback', @(h,e)Fig_preview_spk_thresh_(hFig));

    mh_edit_filter = uimenu(mh_edit, 'Label', 'Filter mode');
    uimenu_options_(mh_edit_filter, {'ndiff', 'bandpass', 'sgdiff', 'fir1', 'user'}, @Fig_preview_filter_, hFig);
    menu_checkbox_(mh_edit_filter, get_filter_(P));

    mh_edit_ref = uimenu(mh_edit, 'Label', 'Reference mode');
    uimenu_options_(mh_edit_ref, {'none', 'mean', 'median'}, @Fig_preview_ref_, hFig); % @TODO: local mean
    menu_checkbox_(mh_edit_ref, P.vcCommonRef);

    uimenu(mh_edit, 'Label', 'Common reference threshold', 'Callback', @(h,e)Fig_preview_ref_thresh_(hFig));

    uimenu(mh_edit, 'Label', 'FFT cleanup threshold', 'Callback', @(h,e)Fig_preview_fft_thresh_(hFig));


    % View menu
    mh_view = uimenu(hFig, 'Label','View');

    mh_view_trange = uimenu(mh_view, 'Label', 'Display time range (s)');
    uimenu_options_(mh_view_trange, {'0.05', '0.1', '0.2', '0.5', '1', '2', '5', 'Custom'}, @Fig_preview_trange_, hFig);
    menu_checkbox_(mh_view_trange, '0.1');
    uimenu(mh_view, 'Label', 'Display site range', 'Callback', @(h,e)Fig_preview_site_range_(hFig));

    uimenu(mh_view, 'Label', 'Show raw traces [F]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'f'), 'Tag', 'menu_preview_view_filter');
    uimenu(mh_view, 'Label', 'Show spike [T]hreshold', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 't'), 'Tag', 'menu_preview_view_threshold');
    uimenu(mh_view, 'Label', 'Hide [S]pikes', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'), 'Tag', 'menu_preview_view_spike');
    uimenu(mh_view, 'Label', 'Hide [G]rid', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'g'), 'Tag', 'menu_preview_view_grid');

    mh_view_site = uimenu(mh_view, 'Label', 'Site view');
    uimenu_options_(mh_view_site, {'Site correlation', 'Spike threshold', 'Event rate (Hz)', 'Event SNR (median)'}, @Fig_preview_site_plot_, hFig);
    menu_checkbox_(mh_view_site, 'Site correlation');

    mh_view_ref = uimenu(mh_view, 'Label', 'Reference view');
    uimenu_options_(mh_view_ref, {'original', 'binned'}, @Fig_preview_view_ref_, hFig);
    menu_checkbox_(mh_view_ref, 'original');

    mh_view_freq = uimenu(mh_view, 'Label', 'Frequency scale');
    uimenu_options_(mh_view_freq, {'Linear', 'Log'}, @Fig_preview_psd_plot_, hFig);
    menu_checkbox_(mh_view_freq, 'Linear');

    % mh_view_psd = uimenu(mh_view, 'Label', 'PSD view');
    % uimenu_options_(mh_view_psd, {'original', 'detrended'}, @Fig_preview_view_psd_, hFig);
    % menu_checkbox_(mh_view_psd, 'original');

    % 'Power', 'Detrended',
end %func
