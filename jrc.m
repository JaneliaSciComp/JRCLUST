% JRCLUST v4
% Alan Liddell, Vidrio Technologies
% Originally written by James Jun
function hJRC_ = jrc(varargin)
    %JRC
    if nargout > 0
        hJRC_ = [];
    end

    if nargin < 1
        disp(jrclust.utils.help()); % later: GUI
        return;
    end

    % discard try-catch here for now as these messages are next to useless
    hJRC = JRC(varargin{:});

    if hJRC.inProgress()
        hJRC.hCfg.verbose = 1; % invocation from command line is always verbose
        hJRC.run();
    end

    if hJRC.isError
        error(hJRC.errMsg);
    end

    if nargout > 0
        hJRC_ = hJRC;
    elseif ~isempty(hJRC.res) % export results to workspace
        jrclust.utils.exportToWorkspace(struct('res', hJRC.res), 0);
    end

    return;
end

%%
%     case 'import-tsf'
%         import_tsf_(vcArg1);
%
%     case 'import-h5'
%         import_h5_(vcArg1);
%
%     case 'import-intan'
%         vcFile_prm_ = import_intan_(vcArg1, vcArg2, vcArg3);
%         return;
%
%     case {'import-nsx', 'import-ns5'}
%         vcFile_prm_ = import_nsx_(vcArg1, vcArg2, vcArg3);
%         return;
%
%     case 'nsx-info'
%         [~, ~, S_file] = nsx_info_(vcArg1);
%         assignWorkspace_(S_file);
%         return;
%
%     case 'load-nsx'
%         load_nsx_(vcArg1);
%         return;
%
%     case 'import-gt'
%         import_gt_silico_(vcArg1);
%
%     case 'unit-test'
%         unit_test_(vcArg1);
%
%     case 'test'
%         varargout{1} = test_(vcArg1, vcArg2, vcArg3, vcArg4, vcArg5);
%
%     case {'make-trial', 'maketrial', 'load-trial', 'loadtrial'}
%         make_trial_(vcFile_prm, 0);
%
%     case {'loadtrial-imec', 'load-trial-imec', 'make-trial-imec', 'maketrial-imec'}
%         make_trial_(vcFile_prm, 1);
%
%     case 'batch'
%         batch_(vcArg1, vcArg2);
%
%     case {'batch-verify', 'batch-validate'}
%         batch_verify_(vcArg1, vcArg2);
%
%     case {'batch-plot', 'batch-activity'}
%         batch_plot_(vcArg1, vcArg2);
%
%     case 'describe'
%         describe_(vcFile_prm);
%
%     case {'import-kilosort', 'import-ksort'}
%         import_ksort_(vcFile_prm);
%
%     case 'import-silico'
%         import_silico_(vcFile_prm, 0);
%
%     case 'import-silico-sort'
%         import_silico_(vcFile_prm, 1);
%
%     case 'export-imec-sync'
%         export_imec_sync_(vcFile_prm);
%
%     case 'export-prm'
%         export_prm_(vcFile_prm, vcArg2);
%
%     case 'preview-test'
%         preview_(P, 1);
%         gui_test_(P, 'Fig_preview');
%
%     case 'traces-lfp'
%         traces_lfp_(P)
%
%     case 'traces-test'
%         traces_(P);
%         traces_test_(P);
%
%     case 'manual-test'
%         manual_(P, 'debug');
%         manual_test_(P);
%         return;
%
%     case 'manual-test-menu'
%         manual_(P, 'debug');
%         manual_test_(P, 'Menu');
%         return;
%
%     case 'export-spk'
%         S0 = get(0, 'UserData');
%         trSpkWav = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkwav.jrc'), 'int16', S0.dimm_spk);
%         assignWorkspace_(trSpkWav);
%
%     case 'export-raw'
%         S0 = get(0, 'UserData');
%         trWav_raw = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkraw.jrc'), 'int16', S0.dimm_spk);
%         assignWorkspace_(trWav_raw);
%
%     case {'export-spkwav', 'spkwav'}
%         export_spkwav_(P, vcArg2); % export spike waveforms
%
%     case {'export-chan'}
%         export_chan_(P, vcArg2); % export channels
%
%     case {'export-car'}
%         export_car_(P, vcArg2); % export common average reference
%
%     case {'export-spkwav-diff', 'spkwav-diff'}
%         export_spkwav_(P, vcArg2, 1); % export spike waveforms
%
%     case 'export-spkamp'
%         export_spkamp_(P, vcArg2); %export microvolt unit
%
%     case {'export-csv', 'exportcsv'}
%         export_csv_(P);
%
%     case {'export-quality', 'exportquality', 'quality'}
%         export_quality_(P);
%
%     case {'export-csv-msort', 'exportcsv-msort'}
%         export_csv_msort_(P);
%
%     case {'export-fet', 'export-features', 'export-feature'}
%         export_fet_(P);
%
%     case 'export-diff'
%         export_diff_(P); %spatial differentiation for two column probe
%
%     case 'import-lfp'
%         import_lfp_(P);
%
%     case 'export-lfp'
%         export_lfp_(P);
%
%     case 'drift'
%         plot_drift_(P);
%
%     case 'plot-rd'
%         plot_rd_(P);
%
%     if contains_(lower(vcCmd), {'verify', 'validate'})
%         if ~is_detected_(P)
%             detect_(P);
%         end
%         if ~is_sorted_(P)
%             sort_(P);
%         end
%         validate_(P);
%     elseif contains_(lower(vcCmd), {'filter'})
%         TWfilter_(P);
%     elseif fError
%         help_();
%     end
