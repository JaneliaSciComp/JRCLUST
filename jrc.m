%--------------------------------------------------------------------------
% JRCLUST v3
% James Jun

function varargout = jrc(cmd, varargin)
    % Memory-efficient version
    % P is static and loaded from file
    % Dynamic variables are set in S0=get(0,'UserData')

    % Add paths
    [dirname, ~] = fileparts(fullfile(mfilename('fullpath')));
    addpath(fullfile(dirname, 'meta')); % info functions
    addpath(fullfile(dirname, 'filesystem')); % file-related functions
    addpath(fullfile(dirname, 'probe')); % probe-related functions
    addpath(fullfile(dirname, 'neuroshare')); % NS-related functions
    addpath(fullfile(dirname, 'struct')); % general struct-related functions
    addpath(fullfile(dirname, 'utils')); % miscellaneous (but useful) tools
    addpath(fullfile(dirname, 'params')); % parameter-related functions
    addpath(fullfile(dirname, 'clusters')); % functions for manipulating the cluster struct
    addpath(fullfile(dirname, 'gui')); % GUI functions
    addpath(fullfile(dirname, 'manual')); % manual curation step functions
    addpath(fullfile(dirname, 'kilosort')); % kilosort-related functions

    % process arguments
    [arg1, arg2, arg3, arg4, arg5] = deal('');
    switch nargin
        case 2
            arg1 = varargin{1};
        case 3
            arg1 = varargin{1}; arg2 = varargin{2};
        case 4
            arg1 = varargin{1}; arg2 = varargin{2}; arg3 = varargin{3};
        case 5
            arg1 = varargin{1}; arg2 = varargin{2}; arg3 = varargin{3};
            arg4 = varargin{4};
        case 6
            arg1 = varargin{1}; arg2 = varargin{2}; arg3 = varargin{3};
            arg4 = varargin{4}; arg5 = varargin{5};
        case 0
            cmd = 'help';
    end

    warning off;

    %-----
    % Command type A: supporting functions
    doExit = 1;

    switch lower(cmd)
        % No arguments
        case 'version'
            jrcVersion(arg1);
        case {'help', '-h', '?', '--help'}
            help_(arg1);
            about_();
        case 'about'
            about_();
        case 'clear'
            clear_(arg1);
        % case 'gui' % TODO: reference 'paramFile' before definition
        %     gui_(arg1, paramFile);
        case {'import-kilosort', 'import-ksort'}
            importKiloSort(arg1, arg2);
        case {'makeprm', 'createprm', 'makeprm-all'}
            paramFile = makeprm_(arg1, arg2, 1, arg3);
            if nargout > 0
                varargout{1} = paramFile;
            end

            if isempty(paramFile)
                return;
            end

            if strcmp(cmd, 'makeprm-all')
                jrc('all', paramFile);
            end
        case 'makeprm-f', makeprm_(arg1, arg2, 0, arg3);
        case 'import-tsf', import_tsf_(arg1);
        case 'import-h5', import_h5_(arg1);
        case 'import-intan', paramFile = import_intan_(arg1, arg2, arg3); return;
        case {'import-nsx', 'import-ns5'}, paramFile = import_nsx_(arg1, arg2, arg3); return;
        case 'nsx-info', [~, ~, S_file] = nsx_info_(arg1); assignWorkspace_(S_file); return;
        case 'load-nsx', load_nsx_(arg1); return;
        case 'load-bin'
            mnWav = load_bin_(arg1, arg2);
            assignWorkspace_(mnWav);
        case 'import-gt', import_gt_silico_(arg1);
        case 'unit-test', unit_test_(arg1, arg2, arg3);
        case 'compile', compile_cuda_(arg1);
        case 'compile-ksort', compile_ksort_();
        case 'test', varargout{1} = test_(arg1, arg2, arg3, arg4, arg5);
        case 'call', varargout{1} = call_(arg1, arg2, arg3);
        case 'export', export_(arg1, arg2, arg3);
        case {'dependencies', 'toolbox', 'toolboxes'}, disp_dependencies_();
        otherwise, doExit = 0;
    end
    if doExit, return; end

    %-----
    % Command type B: Requires .prm file
    if nargin >= 2
        paramFile = getParamFile(arg1);
    else
        paramFile = getParamFile();
        disp(['Working on ', paramFile])
    end

    % stopgap
    twelve(paramFile);

    doExit = 1;
    switch lower(cmd)
        case 'probe', probe_(paramFile);
        case {'make-trial', 'maketrial', 'load-trial', 'loadtrial'}, make_trial_(paramFile, 0);
        case {'loadtrial-imec', 'load-trial-imec', 'make-trial-imec', 'maketrial-imec'}, make_trial_(paramFile, 1);
        case 'edit', edit_(paramFile);
        case 'batch', batch_(arg1, arg2);
        case {'batch-verify', 'batch-validate'}, batch_verify_(arg1, arg2);
        case {'batch-plot', 'batch-activity'}, batch_plot_(arg1, arg2);
        case 'describe', describe_(paramFile);
        case 'import-silico', import_silico_(paramFile, 0);
        case 'import-silico-sort', import_silico_(paramFile, 1);
        case 'export-imec-sync', export_imec_sync_(paramFile);
        case 'export-prm', exportParams(paramFile, arg2);
        case 'dir'
            if any(paramFile == '*') % handle globs
                dir_files_(paramFile, arg2, arg3);
            else
                doExit = 0;
            end
        otherwise
            doExit = 0;
    end
    if doExit, return; end

    %-----
    % Command type C: Requires P structure (loaded from .prm)
    if ~matchFileExt_(paramFile, '.prm'), fprintf(2, 'Must provide .prm file\n'); return ;end
    if ~fileExists(paramFile), fprintf(2, 'File does not exist: %s\n', paramFile); return ;end
    P = loadParams(paramFile);
    if isempty(P), return; end
    fError = 0;
    switch lower(cmd)
        case 'preview', preview_(P);
        case 'preview-test', preview_(P, 1); gui_test_(P, 'Fig_preview');
        case 'traces', traces_(P, 0, arg2);
        case 'traces-lfp', traces_lfp_(P)
        case 'dir', dir_files_(P.multiFilenames);
        case 'traces-test'
            traces_(P, 1); traces_test_(P);
        case {'full', 'all'}
            fprintf('Performing "jrc detect", "jrc sort", "jrc manual" operations.\n');
            detect_(P);
            sort_(P, 0);
            describe_(P.paramFile);
            manual(P);
            return;
        case {'spikesort', 'detectsort', 'detect-sort', 'spikesort-verify', 'spikesort-validate', 'spikesort-manual', 'detectsort-manual'}
            fprintf('Performing "jrc detect", "jrc sort" operations.\n');
            detect_(P);
            sort_(P, 0);
            describe_(P.paramFile);
        case {'detect', 'spikedetect'}
            detect_(P); describe_(P.paramFile);
        case {'sort', 'cluster', 'clust', 'sort-verify', 'sort-validate', 'sort-manual'}
            if ~is_detected_(P)
                detect_(P); sort_(P,0);
            else
                sort_(P);
            end
            describe_(P.paramFile);
        case {'auto', 'auto-verify', 'auto-manual'}
            auto_(P); describe_(P.paramFile);
        case 'manual-test'
            manual(P, 'debug'); manual_test_(P); return;
        case 'manual-test-menu'
            manual(P, 'debug'); manual_test_(P, 'Menu'); return;
        case {'export-wav', 'wav'} % load raw and assign workspace
            mnWav = load_file_(P.vcFile, [], P);
            assignWorkspace_(mnWav);
        case 'export-spk'
            S0 = get(0, 'UserData');
            trSpkWav = load_bin_(strrep(P.paramFile, '.prm', '_waveforms.bin'), 'int16', S0.waveformDims);
            assignWorkspace_(trSpkWav);
        case 'export-raw'
            S0 = get(0, 'UserData');
            trWav_raw = load_bin_(strrep(P.paramFile, '.prm', '_traces.bin'), 'int16', S0.waveformDims);
            assignWorkspace_(trWav_raw);
        case {'export-spkwav', 'spkwav'}, export_spkwav_(P, arg2); % export spike waveforms
        case {'export-chan'}, export_chan_(P, arg2); % export channels
        case {'export-car'}, export_car_(P, arg2); % export common average reference
        case {'export-spkwav-diff', 'spkwav-diff'}, export_spkwav_(P, arg2, 1); % export spike waveforms
        case 'export-spkamp', export_spkamp_(P, arg2); %export microvolt unit
        case {'export-csv', 'exportcsv'}, export_csv_(P);
        case {'export-quality', 'exportquality', 'quality'}, export_quality_(P);
        case {'export-csv-msort', 'exportcsv-msort'}, export_csv_msort_(P);
        case {'activity', 'plot-activity'}, plot_activity_(P);
        case {'export-fet', 'export-features', 'export-feature'}, export_fet_(P);
        case 'export-diff', export_diff_(P); %spatial differentiation for two column probe
        case 'import-lfp', import_lfp_(P);
        case 'export-lfp', export_lfp_(P);
        case 'drift', plot_drift_(P);
        case 'plot-rd', plot_rd_(P);
        otherwise, fError = 1;
    end %switch

    % supports compound commands (ie. 'sort-verify', 'sort-manual').
    if contains_(lower(cmd), {'verify', 'validate'})
        if ~is_detected_(P), detect_(P); end
        if ~isSorted(P)
            sort_(P);
        end
        validate_(P);
    elseif contains_(lower(cmd), {'manual',' gui', 'ui'})
        manual(P);
    elseif contains_(lower(cmd), {'filter'})
        TWfilter_(P);
    elseif fError
        help_();
    end
end % function
