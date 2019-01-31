classdef JRC < handle & dynamicprops
    %JRC Handle class for spike sorting pipeline

    properties (SetObservable, SetAccess=private, Hidden)
        errMsg;         %
        isCompleted;    %
        isError;        %
        isDetect;       %
        isSort;         %
        isCurate;       %
    end

    properties (SetObservable, SetAccess=private)
        args;           %
        cmd;            %
        hCfg;           %
        hDetect;        %
        hSort;          %
        hCurate;        %
    end

    % most general data from detection step
    properties (Dependent)
        spikeTimes;     %
        spikeSites;     % 
    end

    % data from clustering step
    properties (Dependent)
        hClust;
    end

    % structs containing results from steps in the pipeline
    properties (Hidden, SetAccess=private, SetObservable)
        res;            % results struct
    end

    %% LIFECYCLE
    methods
        function obj = JRC(varargin)
            %JRC Construct an instance of this class
            obj.args = varargin;
            obj.isCompleted = 0;
            obj.isError = 0;
            obj.isDetect = 0;
            obj.isSort = 0;
            obj.isCurate = 0;

            if ~jrclust.utils.sysCheck()
                obj.errMsg = 'system requirements not met';
                obj.isError = 1;
            elseif nargin > 0
                % handle arguments (legacy mode)
                obj.processArgs();
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected)
        function processArgs(obj)
            %PROCESSARGS Handle command-line arguments, load param file
            nargs = numel(obj.args);

            if nargs == 0
                return;
            end

            obj.cmd = lower(obj.args{1});
            obj.args = obj.args(2:end);
            nargs = nargs - 1;

            % paired manual commands
            if contains(obj.cmd, '-manual')
                obj.cmd = strrep(obj.cmd, '-manual', '');
                obj.isCurate = 1;
            end

            switch obj.cmd
                % deprecated commands; will be removed in a future release
                case {'compile-ksort', 'dir', 'edit', 'git-pull', 'issue', 'import-kilosort-sort', ...
                      'import-ksort-sort', 'kilosort', 'kilosort-verify', 'ksort', 'ksort-verify' ...
                      'which', 'wiki', 'wiki-download'}
                    jrclust.utils.depWarn(obj.cmd);
                    obj.isCompleted = 1;
                    return;

                case {'doc', 'doc-edit'}
                    imsg = 'Please visit the wiki at https://github.com/JaneliaSciComp/JRCLUST/wiki';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case 'download'
                    imsg = 'You can find sample.bin and sample.meta at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case 'gui'
                    imsg = 'GUI is not implemented yet';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case 'install'
                    imsg = 'You might be looking for `compile` instead';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case {'set', 'setprm', 'set-prm'}
                    imsg = 'Create a new JRC handle instead';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case 'update'
                    imsg = 'Please check the repository at https://github.com/JaneliaSciComp/JRCLUST for updates';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                case {'load-bin', 'export-wav', 'wav'}
                    imsg = 'Please use jrclust.models.recordings.Recording instead';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = 1;
                    return;

                % deprecated synonyms, warn but proceed
                case 'spikedetect'
                    obj.cmd = 'detect';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);

                case {'cluster', 'clust', 'sort-verify', 'sort-validate'}
                    obj.cmd = 'sort';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);

                case {'detectsort', 'spikesort', 'spikesort-verify', 'spikesort-validate'}
                    obj.cmd = 'detect-sort';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);

                case 'all'
                    obj.cmd = 'full';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);
                
                case 'plot-activity'
                    obj.cmd = 'activity';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);

                case {'makeprm', 'createprm'}
                    obj.cmd = 'bootstrap';
                    imsg = sprintf('Please use ''%s'' in the future', obj.cmd);
                    jrclust.utils.depWarn(obj.cmd, imsg);

                % info commands
                case 'about'
                    md = jrclust.utils.info();
                    verstr = sprintf('%s v%s', md.program, jrclust.utils.version());
                    abstr = jrclust.utils.about();
                    msgbox(abstr, verstr);
                    obj.isCompleted = 1;
                    return;

                case 'help'
                    disp(jrclust.utils.help());
                    obj.isCompleted = 1;
                    return;

                case 'version'
                    md = jrclust.utils.info();
                    fprintf('%s v%s\n', md.program, jrclust.utils.version());
                    obj.isCompleted = 1;
                    return;

                % workflow commands
                case 'bootstrap'
                    doBootstrap(obj.args{:});
                    obj.isCompleted = 1;
                    return;

                % preview commands
                case 'probe'
                    if nargs == 0
                        obj.errMsg = 'Specify a probe file or config file';
                        obj.isError = 1;
                        return;
                    end

                    probeFile = obj.args{1};
                    
                    if endsWith(probeFile, '.prm')
                        hCfg_ = jrclust.Config(probeFile);
                        doPlotProbe(hCfg_);
                    else % not a config file
                        [~, ~, ext] = fileparts(probeFile);
                        if isempty(ext) % a convenience for a forgetful mind
                            probeFile = [probeFile '.prb'];
                        end

                        probeFile_ = jrclust.utils.absPath(probeFile, fullfile(jrclust.utils.basedir(), 'probes'));
                        if isempty(probeFile_)
                            obj.errMsg = sprintf('Could not find probe file: %s', probeFile);
                            obj.isError = 1;
                            return;
                        end

                        probeData = doLoadProbe(probeFile_);
                        doPlotProbe(probeData);
                    end

                    obj.isCompleted = 1;
                    return;

                case 'preview'
                    hCfg_ = jrclust.Config(obj.args{1});
                    hPreview = jrclust.controllers.curate.PreviewController(hCfg_);
                    hPreview.preview();

                    obj.isCompleted = 1;
                    return;

                    
                case 'traces'
                    hCfg_ = jrclust.Config(obj.args{1});
                    hTraces = jrclust.controllers.curate.TracesController(hCfg_);
                    if numel(obj.args) > 1
                        recID = str2double(obj.args{2});
                        if isnan(recID)
                            recID = [];
                        end
                    else
                        recID = [];
                    end
                    hTraces.show(recID, 0, obj.hClust);

                    obj.isCompleted = 1;
                    return;                    
            end

            % command sentinel
            detectCmds = {'detect', 'detect-sort', 'full'};
            sortCmds   = {'sort', 'detect-sort', 'full'};
            curateCmds = {'manual', 'full'};
            miscCmds = {'activity', 'auto'};

            legalCmds = unique([detectCmds, sortCmds curateCmds miscCmds]);

            if ~any(strcmpi(obj.cmd, legalCmds))
                obj.errMsg = sprintf('Command `%s` not recognized', obj.cmd);
                errordlg(obj.errMsg, 'Unrecognized command');
                obj.isError = 1;
                return;
            end

            % determine which commands in the pipeline to run
            if any(strcmp(obj.cmd, curateCmds))
                obj.isCurate = 1;
            end

            if any(strcmp(obj.cmd, sortCmds))
                obj.isSort = 1;
            end

            if any(strcmp(obj.cmd, detectCmds))
                obj.isDetect = 1;
            end

            % commands from here on out require a parameter file
            if nargs < 1
                obj.errMsg = sprintf('Command `%s` requires a parameter file', obj.cmd);
                errordlg(obj.errMsg, 'Missing parameter file');
                obj.isError = 1;
                return;
            end

            % load parameter file
            configFile = obj.args{1};
            obj.hCfg = jrclust.Config(configFile);
        end
    end

    %% USER METHODS
    methods
        function ip = inProgress(obj)
            %INPROGRESS Return true if in progress (not finished or errored)
            ip = ~(obj.isCompleted || obj.isError);
        end

        function it = invocationType(obj)
            if obj.isError
                it = 'error';
            elseif obj.isSort
                it = 'sort';
            elseif obj.isDetect
                it = 'detect';
            elseif obj.isCurate
                it = 'curate';
            else % ~(obj.isCurate || obj.isSort)
                it = 'info';
            end
        end

        function rerun(obj)
            %RERUN Rerun commands
            if obj.isError
                error(obj.errMsg);
            else
                obj.isCompleted = 0;
                obj.run();
            end
        end

        function run(obj)
            %RUN Run commands
            if obj.isError
                error(obj.errMsg);
            elseif obj.isCompleted
                warning('command ''%s'' completed successfully; to rerun, use rerun()', obj.cmd);
                return;
            end

            % try to warm up the local parallel pool before taking a swim
            if obj.hCfg.useParfor && (obj.isDetect || obj.isSort || obj.isSort)
                try
                    parpool('local');
                catch
                end
            end

            % load saved data only if we're not starting over
            if ~obj.isDetect
                obj.res = obj.loadFiles();
            else
                obj.res = [];
            end

            if obj.isCurate && ~obj.isSort
                if ~isfield(obj.res, 'hClust')
                    dlgAns = questdlg('Could not find all required data. Sort?', 'Sorting required', 'No');
                    if strcmp(dlgAns, 'Yes')
                        obj.isSort = 1;
                    else
                        return;
                    end
                end
            end

            if obj.isSort && ~obj.isDetect
                if ~isfield(obj.res, 'spikeTimes')
                    obj.isDetect = 1;
                end
            end

            doSave = obj.isSort || obj.isDetect;

            % notify user that we're using previously-computed results
            if obj.hCfg.verbose
                if ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
                    fprintf('Using spikes detected on %s\n', datestr(obj.res.detectedOn));
                elseif ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
                    fprintf('Using previously-detected spikes\n');
                end
                if ~isempty(obj.res) && isfield(obj.res, 'sortedOn')
                    fprintf('Using clustering computed on %s\n', datestr(obj.res.sortedOn));
                elseif ~isempty(obj.res) && isfield(obj.res, 'hClust')
                    fprintf('Using previously-clustered spikes\n');
                end
                if ~isempty(obj.res) && isfield(obj.res, 'curatedOn')
                    fprintf('Last manually edit on %s\n', datestr(obj.res.curatedOn));
                end
            end

            % set random seeds
            rng(obj.hCfg.randomSeed);
            if obj.hCfg.useGPU
                % while we're here, clear GPU memory
                if obj.isDetect || obj.isSort
                    if obj.hCfg.verbose
                        fprintf('Clearing GPU memory...');
                    end
                    gpuDevice(); % selects GPU device
                    gpuDevice([]); % clears GPU memory
                    if obj.hCfg.verbose
                        fprintf('done\n');
                    end
                end

                parallel.gpu.rng(obj.hCfg.randomSeed);
            end

            % PLOT ACTIVITY OR RECLUSTER
            if strcmp(obj.cmd, 'activity')
                if ~isempty(obj.res) && all(cellfun(@(f) ismember(f, fieldnames(obj.res)), ...
                                                    {'spikeTimes', 'spikeSites', 'spikeAmps'}))
                    doPlotActivity(obj.hCfg, obj.res);
                end
            elseif strcmp(obj.cmd, 'auto')
                if ~isempty(obj.res) && isfield(obj.res, 'hClust')
                    obj.res.hClust.hCfg = obj.hCfg; % update hClust's config
                    obj.res.hClust.reassign();
                    obj.res.hClust.autoMerge();
                    obj.res.sortedOn = now();
                    doSave = 1;
                else
                    obj.isError = 1;
                    obj.errMsg = 'hClust not found';
                    return;
                end
            end

            % DETECT SPIKES
            gpuDetect = obj.hCfg.useGPU; % save this in case useGPU is disabled during detection step
            if obj.isDetect
                obj.hDetect = jrclust.controllers.detect.DetectController(obj.hCfg);
                dRes = obj.hDetect.detect();

                if obj.hDetect.isError
                    error(obj.hDetect.errMsg);
                elseif obj.hCfg.verbose
                    fprintf('Detection completed in %0.2f seconds\n', dRes.runtime);
                end

                obj.res = dRes;
            end

            % CLUSTER SPIKES
            gpuSort = obj.hCfg.useGPU | gpuDetect;
            if obj.isSort
                obj.hCfg.useGPU = gpuSort;

                obj.hSort = jrclust.controllers.sort.SortController(obj.hCfg);
                sRes = obj.hSort.sort(obj.res);

                if obj.hSort.isError
                    error(obj.hSort.errMsg);
                elseif obj.hCfg.verbose
                    fprintf('Sorting completed in %0.2f seconds\n', sRes.runtime);
                end

                obj.res = jrclust.utils.mergeStructs(obj.res, sRes);
            end

            % save our results for later
            if doSave
                obj.saveFiles();
            end

            % CURATE SPIKES
            gpuCurate = obj.hCfg.useGPU | gpuSort;
            if obj.isCurate
                obj.hCfg.useGPU = gpuCurate;

                obj.hCurate = jrclust.controllers.curate.CurateController(obj.hClust);
                obj.hCurate.beginSession();
            end

            obj.isCompleted = 1;
        end

        function saveFiles(obj)
            %SAVEFILES Save results struct to disk
            if obj.isError
                error(obj.errMsg);
            end

            doSaveFiles(obj.res, obj.hCfg);
        end

        function res = loadFiles(obj)
            %LOADFILES Load results struct
            if obj.isError
                error(obj.errMsg);
            end

            res = doLoadFiles(obj.hCfg);
        end
    end

    % GETTERS/SETTERS
    methods
        % hClust
        function hc = get.hClust(obj)
            if isempty(obj.res) || ~isfield(obj.res, 'hClust')
                hc = [];
            else
                hc = obj.res.hClust;
            end
        end

        % isError
        function ie = get.isError(obj)
            ie = obj.isError;
            if ~isempty(obj.hCfg) % pick up error in Config
                ie = ie || obj.hCfg.isError;
            end
            if ~isempty(obj.hDetect) % pick up error in hDetect
                ie = ie || obj.hDetect.isError;
            end
            if ~isempty(obj.hSort) % pick up error in hSort
                ie = ie || obj.hSort.isError;
            end
        end

        % spikeTimes
        function st = get.spikeTimes(obj)
            if isempty(obj.res)
                st = [];
            else
                st = obj.res.spikeTimes;
            end
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            if isempty(obj.res)
                ss = [];
            else
                ss = obj.res.spikeSites;
            end
        end
    end
end

