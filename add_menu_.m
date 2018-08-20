%--------------------------------------------------------------------------
function add_menu_(hFig, P)
    drawnow;
    posvec = get(hFig, 'OuterPosition');

    set(hFig, 'MenuBar','None');
    mh_file = uimenu(hFig,'Label','File');
    uimenu(mh_file,'Label', 'Save', 'Callback', @save_manual_);
    uimenu(mh_file,'Label', 'Save figures as .fig', 'Callback', @(h,e)save_figures_('.fig'));
    uimenu(mh_file,'Label', 'Save figures as .png', 'Callback', @(h,e)save_figures_('.png'));
    uimenu(mh_file,'Label', 'Describe', 'Callback', @(h,e)msgbox_(describe_()), 'Separator', 'on');
    uimenu(mh_file,'Label', 'Edit prm file', 'Callback', @edit_prm_);
    uimenu(mh_file,'Label', 'Reload prm file', 'Callback', @reload_prm_);
    uimenu(mh_file,'Label', 'Export units to csv', 'Callback', @export_csv_, 'Separator', 'on');
    uimenu(mh_file,'Label', 'Export unit qualities to csv', 'Callback', @(h,e)export_quality_);
    uimenu(mh_file,'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
    uimenu(mh_file,'Label', 'Export selected mean unit waveforms', 'Callback', @(h,e)export_mrWav_clu_);
    uimenu(mh_file,'Label', 'Export all waveforms from the selected unit', 'Callback', @(h,e)export_spikeWaveforms_);
    uimenu(mh_file,'Label', 'Export firing rate for all units', 'Callback', @(h,e)export_rate_);
    uimenu(mh_file,'Label', 'Exit', 'Callback', @exit_manual_, 'Separator', 'on', 'Accelerator', 'Q');

    mh_edit = uimenu(hFig,'Label','Edit');
    uimenu(mh_edit,'Label', '[M]erge', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'm'));
    uimenu(mh_edit,'Label', 'Merge auto', 'Callback', @(h,e)merge_auto_());
    uimenu(mh_edit,'Label', '[D]elete', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'd'), 'Separator', 'on');
    uimenu(mh_edit,'Label', 'Delete auto', 'Callback', @(h,e)delete_auto_());
    uimenu(mh_edit,'Label', 'Delete annotated', 'Callback', @(h,e)delete_annotate()); % TW
    uimenu(mh_edit,'Label', '[S]plit', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'), 'Separator', 'on');
    uimenu(mh_edit,'Label', 'Auto split max-chan', 'Callback', @(h,e)auto_split_(0));
    uimenu(mh_edit,'Label', 'Auto split multi-chan', 'Callback', @(h,e)auto_split_(1));
    uimenu(mh_edit,'Label', 'Annotate', 'Callback', @(h,e)unit_annotate_());

    mh_view = uimenu(hFig,'Label','View');
    uimenu(mh_view,'Label', 'Show traces', 'Callback', @(h,e)traces_());
    uimenu(mh_view,'Label', 'View all [R]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'r'));
    uimenu(mh_view,'Label', '[Z]oom selected', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'z'));
    uimenu(mh_view,'Label', '[W]aveform (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'w'));
    uimenu(mh_view,'Label', '[N]umbers (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'n'));
    uimenu(mh_view,'Label', 'Show raw waveform', 'Callback', @(h,e) showRawWaveforms(h), ...
        'Checked', ifeq_(get_(P, 'fWav_raw_show'), 'on', 'off'));

    %uimenu(mh_view,'Label', 'Threshold by sites', 'Callback', @(h,e)keyPressFcn_thresh_(hFig, 'n'));
    % uimenu(mh_view,'Label', '.prm file', 'Callback', @edit_prm_);
    uimenu(mh_view,'Label', 'Reset window positions', 'Callback', @reset_position_);

    mh_proj = uimenu(hFig,'Label','Projection');
    uimenu(mh_proj, 'Label', 'vpp', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.displayFeature, {'vpp', 'vmin'}));
    uimenu(mh_proj, 'Label', 'pca', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.displayFeature, {'pca'}));
    uimenu(mh_proj, 'Label', 'ppca', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.displayFeature, {'ppca', 'private pca'}));
    % uimenu(mh_proj, 'Label', 'cov', 'Callback', @(h,e)proj_view_(h), ...
    %     'Checked', if_on_off_(P.displayFeature, {'cov', 'spacetime'}));

    mh_plot = uimenu(hFig,'Label','Plot');
    uimenu(mh_plot, 'Label', 'All unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_);
    uimenu(mh_plot, 'Label', 'Selected unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_(1));

    mh_info = uimenu(hFig,'Label','','Tag', 'mh_info');
    uimenu(mh_info, 'Label', 'Annotate unit', 'Callback', @unit_annotate_);
    uimenu(mh_info, 'Label', 'single', 'Callback', @(h,e)unit_annotate_(h,e,'single'));
    uimenu(mh_info, 'Label', 'multi', 'Callback', @(h,e)unit_annotate_(h,e,'multi'));
    uimenu(mh_info, 'Label', 'noise', 'Callback', @(h,e)unit_annotate_(h,e,'noise'));
    uimenu(mh_info, 'Label', 'clear annotation', 'Callback', @(h,e)unit_annotate_(h,e,''));
    uimenu(mh_info, 'Label', 'equal to', 'Callback', @(h,e)unit_annotate_(h,e,'=%d'));

    mh_history = uimenu(hFig, 'Label', 'History', 'Tag', 'mh_history');

    mh_help = uimenu(hFig,'Label','Help');
    uimenu(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);

    drawnow;
    set(hFig, 'OuterPosition', posvec);
end %func
