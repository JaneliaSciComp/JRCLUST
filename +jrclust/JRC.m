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
        args;           % arguments passed after command
        cmd;            % command passed
        hDetect;        %
        hSort;          %
        hCurate;        %
    end

    properties (Hidden, SetObservable)
        hCfg;           % Config object
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
                obj.errMsg = 'System requirements not met';
                obj.isError = 1;
            elseif nargin > 0
                try
                    obj.hCfg = varargin{1};
                catch ME % arg 1 is not a Config, process arguments
                    % handle arguments (legacy mode)
                    obj.processArgs();
                end
            end
        end
    end

    %% UTILITY METHODS
    methods (Access=protected)
        function clearMemory(obj)
            %CLEARMEMORY Clear GPU memory and set random seeds
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

            rng(obj.hCfg.randomSeed);
        end

        function deprecateCmd(obj, oldCmd, iMsg, newCmd)
            %DEPRECATECMD Print a message regarding deprecation
            %   if newCmd is specified, sets obj.cmd to newCmd
            if nargin < 3
                iMsg = '';
            end
            if nargin < 4
                newCmd = '';
            end

            % set a default iMsg in case of an alias
            if isempty(iMsg) && ~isempty(newCmd)
                iMsg = sprintf('Please use ''%s'' in the future', newCmd);
            end

            jrclust.utils.depWarn(oldCmd, iMsg);

            % set obj.cmd here
            if ~isempty(newCmd)
                obj.cmd = newCmd;
            end
        end

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

            if nargs > 0
                % load parameter file
                configFile = obj.args{1};
                try
                    obj.hCfg = jrclust.Config(configFile);
                    % save imported config file
                    if obj.hCfg.isV3Import
                        obj.hCfg.save();
                    end
                catch ME % arg 1 need not be a config file
                end
            end

            switch obj.cmd
                % deprecated commands; will be removed in a future release
                case {'compile-ksort', 'dir', 'edit', 'git-pull', 'issue', 'import-kilosort-sort', ...
                      'import-ksort-sort', 'kilosort', 'kilosort-verify', 'ksort', 'ksort-verify' ...
                      'which', 'wiki', 'wiki-download'}
                    obj.deprecateCmd(obj.cmd);
                    obj.isCompleted = 1;

                case {'doc', 'doc-edit'}
                    iMsg = 'Please visit the wiki at https://github.com/JaneliaSciComp/JRCLUST/wiki';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case 'download'
                    iMsg = 'You can find sample.bin and sample.meta at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case 'gui'
                    iMsg = 'GUI is not implemented yet';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case 'install'
                    iMsg = 'You might be looking for `compile` instead';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case {'set', 'setprm', 'set-prm'}
                    iMsg = 'Create a new JRC handle instead';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case 'update'
                    iMsg = 'Please check the repository at https://github.com/JaneliaSciComp/JRCLUST for updates';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                case {'load-bin', 'export-wav', 'wav'}
                    iMsg = 'Please use jrclust.models.recordings.Recording instead';
                    obj.deprecateCmd(obj.cmd, iMsg);
                    obj.isCompleted = 1;

                % deprecated synonyms, warn but proceed
                case 'spikedetect'
                    obj.deprecateCmd(obj.cmd, '', 'detect');

                case {'cluster', 'clust', 'sort-verify', 'sort-validate'}
                    obj.deprecateCmd(obj.cmd, '', 'sort');

                case {'detectsort', 'spikesort', 'spikesort-verify', 'spikesort-validate'}
                    obj.deprecateCmd(obj.cmd, '', 'detect-sort');

                case 'all'
                    obj.deprecateCmd(obj.cmd, '', 'full');

                case 'auto'
                    obj.deprecateCmd(obj.cmd, '', 'recluster');
                    obj.recluster();
                    obj.isCompleted = ~obj.isError;

                case 'plot-activity'
                    obj.deprecateCmd(obj.cmd, '', 'activity');
                    obj.activity();
                    obj.isCompleted = 1;

                case {'makeprm', 'createprm'}
                    obj.deprecateCmd(obj.cmd, '', 'bootstrap');
                    doBootstrap(obj.args{:});
                    obj.isCompleted = 1;

                % info commands
                case 'about'
                    md = jrclust.utils.info();
                    verstr = sprintf('%s v%s', md.program, jrclust.utils.version());
                    abstr = jrclust.utils.about();
                    msgbox(abstr, verstr);
                    obj.isCompleted = 1;

                case 'help'
                    jrclust.utils.help();
                    obj.isCompleted = 1;

                case 'version'
                    md = jrclust.utils.info();
                    fprintf('%s v%s\n', md.program, jrclust.utils.version());
                    obj.isCompleted = 1;

                % workflow commands
                case 'bootstrap'
                    obj.hCfg = doBootstrap(obj.args{:});
                    obj.isCompleted = 1;

                case 'compile'
                    jrclust.CUDA.compileCUDA();
                    obj.isCompleted = 1;

                case 'importv3'
                    jrclust.import.importv3(obj.args{1});
                    obj.isCompleted = 1;

                % misc commands
                case 'activity'
                    obj.activity();
                    obj.isCompleted = 1;

                case 'probe'
                    if isempty(obj.hCfg) && nargs == 0
                        obj.isError = 1;
                        obj.errMsg = 'Specify a probe or a config file';
                    elseif isempty(obj.hCfg)
                        obj.probe(obj.args{1});
                        obj.isCompleted = 1;
                    else
                        obj.probe();
                        obj.isCompleted = 1;
                    end

                case 'preview'
                    obj.preview();
                    obj.isCompleted = 1;

                case 'recluster'
                    obj.recluster();
                    obj.isCompleted = ~obj.isError;

                case 'traces'
                    obj.traces();
                    obj.isCompleted = 1;
            end

            if ~obj.inProgress
                return;
            end

            % command sentinel
            detectCmds = {'detect', 'detect-sort', 'full'};
            sortCmds   = {'sort', 'detect-sort', 'full'};
            curateCmds = {'manual', 'full'};

            legalCmds = unique([detectCmds, sortCmds curateCmds]);

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
        end

        function startParPool(obj)
            %STARTPARPOOL Start the parallel pool
            if obj.hCfg.useParfor
                try
                    parpool('local');
                catch % parpool already running
                end
            end
        end
    end

    %% PREVIEW/RECLUSTER METHODS
    methods
        function activity(obj)
            %ACTIVITY Plot activity
            if isempty(obj.res)
                obj.loadFiles();
            end

            doPlotActivity(obj.res, obj.hCfg);
        end

        function preview(obj)
            %PREVIEW Display Preview GUI
            hPreview = jrclust.controllers.curate.PreviewController(obj.hCfg);
            hPreview.preview();
        end

        function probe(obj, probeFile)
            %PROBE Plot a probe layout
            if nargin > 1
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
            elseif isempty(obj.hCfg) 
                obj.errMsg = 'Specify a probe file or config file';
                obj.isError = 1;
                return;
            else
                doPlotProbe(obj.hCfg);
            end
        end

        function traces(obj)
            %TRACES Show traces
            if isempty(obj.res)
                obj.loadFiles();
            end

            hTraces = jrclust.controllers.curate.TracesController(obj.hCfg);
            if numel(obj.args) > 1
                recID = str2double(obj.args{2});
                if isnan(recID)
                    recID = [];
                end
            else
                recID = [];
            end

            hTraces.show(recID, 0, obj.hClust);
        end
    end

    %% PIPELINE METHODS
    methods
        function recluster(obj)
            %RECLUSTER Recluster spikes
            if isempty(obj.res)
                obj.loadFiles();
            end

            if isfield(obj.res, 'hClust')
                obj.res.hClust.hCfg = obj.hCfg; % update hClust's config
                obj.res.hClust.reassign();
                obj.res.hClust.autoMerge();
                obj.res.sortedOn = now();
                obj.saveRes(1);
            else
                obj.isError = 1;
                obj.errMsg = 'hClust not found';
            end
        end

        function dRes = detect(obj)
            %DETECT Detect spikes in recording
            obj.res = []; % reset obj.res
            obj.isDetect = 1;

            % clear GPU memory and set random seeds
            obj.clearMemory();

            % start the parallel pool
            if obj.hCfg.useParfor
                obj.startParPool();
            end

            obj.hDetect = jrclust.controllers.detect.DetectController(obj.hCfg);
            dRes = obj.hDetect.detect();

            if obj.hDetect.isError
                error(obj.hDetect.errMsg);
            elseif obj.hCfg.verbose
                fprintf('Detection completed in %0.2f seconds\n', dRes.detectTime);
            end

            obj.res = dRes;

            % save files
            obj.saveBinaries();
            obj.saveRes(0);
        end

        function sRes = sort(obj)
            %SORT Cluster detected spikes
            if isempty(obj.res)
                obj.loadFiles();
            end
            obj.isSort = 1;

            if ~isfield(obj.res, 'spikeFeatures')
                dlgAns = questdlg('Could not find all required data. Detect?', 'Detection required', 'No');
                if strcmp(dlgAns, 'Yes')
                    obj.detect();
                else
                    obj.isCompleted = 1;
                    return;
                end
            end

            if obj.hCfg.verbose
                % inform user we're using previously detected spikes
                if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
                    fprintf('Using spikes detected on %s\n', datestr(obj.res.detectedOn));
                elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
                    fprintf('Using previously-detected spikes\n');
                end
            end

            % clear GPU memory and set random seeds
            obj.clearMemory();

            % start the parallel pool
            if obj.hCfg.useParfor
                obj.startParPool();
            end

            obj.hSort = jrclust.controllers.sort.SortController(obj.hCfg);
            sRes = obj.hSort.sort(obj.res);

            if obj.hSort.isError
                error(obj.hSort.errMsg);
            elseif obj.hCfg.verbose
                fprintf('Sorting completed in %0.2f seconds\n', sRes.sortTime);
            end

            obj.res = jrclust.utils.mergeStructs(obj.res, sRes);
            obj.saveRes(obj.isDetect); % force overwrite if we're detecting
        end

        function curate(obj)
            %CURATE Spin up the manual GUI for curation
            if isempty(obj.res)
                obj.loadFiles();
            end
            obj.isCurate = 1;

            if ~isfield(obj.res, 'hClust')
                dlgAns = questdlg('Could not find all required data. Sort?', 'Sorting required', 'No');
                if strcmp(dlgAns, 'Yes')
                    obj.sort();
                else
                    obj.isCompleted = 1;
                    return;
                end
            end

            if obj.hCfg.verbose
                % inform user we're using previously detected spikes
                if ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'detectedOn')
                    fprintf('Using spikes detected on %s\n', datestr(obj.res.detectedOn));
                elseif ~obj.isDetect && ~isempty(obj.res) && isfield(obj.res, 'spikeTimes')
                    fprintf('Using previously-detected spikes\n');
                end

                % inform user we're using a previously-computed clustering
                if ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'sortedOn')
                    fprintf('Using clustering computed on %s\n', datestr(obj.res.sortedOn));
                elseif ~obj.isSort && ~isempty(obj.res) && isfield(obj.res, 'hClust')
                    fprintf('Using previously-clustered spikes\n');
                end

                % inform user of the last time this set was curated
                if ~isempty(obj.res) && isfield(obj.res, 'curatedOn')
                    fprintf('Last manually edited on %s\n', datestr(obj.res.curatedOn));
                end
            end

            % clear GPU memory and set random seeds
            obj.clearMemory();

            % start the parallel pool
            if obj.hCfg.useParfor
                obj.startParPool();
            end

            obj.hCurate = jrclust.controllers.curate.CurateController(obj.hClust);
            obj.hCurate.beginSession();
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

            % DETECT SPIKES
            gpuDetect = obj.hCfg.useGPU; % save this in case useGPU is disabled during detection step
            if obj.isDetect
                obj.detect();
            end

            % CLUSTER SPIKES
            gpuSort = obj.hCfg.useGPU | gpuDetect;
            if obj.isSort
                obj.hCfg.useGPU = gpuSort;
                obj.sort();
            end

            % CURATE SPIKES
            gpuCurate = obj.hCfg.useGPU | gpuSort;
            if obj.isCurate
                obj.hCfg.useGPU = gpuCurate;
                obj.curate();
            end

            obj.isCompleted = 1;
        end

        function saveBinaries(obj)
            %SAVEBINARIES Save raw/filtered traces, features to disk
            if isempty(obj.res)
                return;
            end
            sessionName = obj.hCfg.sessionName;

            % save spikesRaw
            if isfield(obj.res, 'spikesRaw') && ~isempty(obj.res.spikesRaw)
                filename = fullfile(obj.hCfg.outputDir, [sessionName, '_raw.jrc']);
                [fid, fErr] = fopen(filename, 'w');
                if fid == -1
                    warning('Failed to save spikesRaw: %s', fErr);
                else
                    ct = fwrite(fid, obj.res.spikesRaw, '*int16');
                    fclose(fid);

                    if ct == numel(obj.res.spikesRaw) && obj.hCfg.verbose
                        fprintf('Saved spikesRaw to %s\n', filename);
                    end
                end
            end

            % save spikesFilt
            if isfield(obj.res, 'spikesFilt') && ~isempty(obj.res.spikesFilt)
                filename = fullfile(obj.hCfg.outputDir, [sessionName, '_filt.jrc']);
                [fid, fErr] = fopen(filename, 'w');
                if fid == -1
                    warning('Failed to save spikesFilt: %s', fErr);
                else
                    ct = fwrite(fid, obj.res.spikesFilt, '*int16');
                    fclose(fid);

                    if ct == numel(obj.res.spikesFilt) && obj.hCfg.verbose
                        fprintf('Saved spikesFilt to %s\n', filename);
                    end
                end
            end

            % save spikeFeatures
            if isfield(obj.res, 'spikeFeatures') && ~isempty(obj.res.spikeFeatures)
                filename = fullfile(obj.hCfg.outputDir, [sessionName, '_features.jrc']);
                [fid, fErr] = fopen(filename, 'w');
                if fid == -1
                    warning('Failed to save spikeFeatures: %s', fErr);
                else
                    ct = fwrite(fid, obj.res.spikeFeatures, '*single');
                    fclose(fid);

                    if ct == numel(obj.res.spikeFeatures) && obj.hCfg.verbose
                        fprintf('Saved spikeFeatures to %s\n', filename);
                    end
                end
            end
        end

        function saveRes(obj, forceOverwrite)
            %SAVERES Save results struct
            if isempty(obj.res)
                return;
            end
            if nargin < 2
                forceOverwrite = 0;
            end
            forceOverwrite = forceOverwrite | obj.hCfg.getOr('testRun', 0);

            sessionName = obj.hCfg.sessionName;
            filename = fullfile(obj.hCfg.outputDir, [sessionName '_res.mat']);
             % don't overwrite unless explicitly permitted (or testing)
            if exist(filename, 'file') && ~forceOverwrite
                question = sprintf('%s already exists. Overwrite?', filename);
                dlgAns = questdlg(question, 'Confirm overwrite', 'No');
                if isempty(dlgAns) || ismember(dlgAns, {'No', 'Cancel'})
                    return;
                end
            end

            % save everything else (don't save spikesRaw, spikesFilt,
            % spikeFeatures inside hClust)
            if isfield(obj.res, 'hClust') && ~isempty(obj.res.hClust.spikesRaw)
                doRestore = 1;

                spikesRaw = obj.res.hClust.spikesRaw;
                obj.res.hClust.spikesRaw = [];

                spikesFilt = obj.res.hClust.spikesFilt;
                obj.res.hClust.spikesFilt = [];

                spikeFeatures = obj.res.hClust.spikeFeatures;
                obj.res.hClust.spikeFeatures = [];
            else
                doRestore = 0;
            end

            jrclust.utils.saveStruct(obj.res, filename);

            if doRestore % restore spikesRaw, spikesFilt, spikeFeatures to hClust
                obj.res.hClust.spikesRaw = spikesRaw;
                obj.res.hClust.spikesFilt = spikesFilt;
                obj.res.hClust.spikeFeatures = spikeFeatures;
            end

            if obj.hCfg.verbose
                fprintf('Saved results to %s\n', filename);
            end
        end

        function loadFiles(obj)
            %LOADFILES Load results struct
            if obj.isError
                error(obj.errMsg);
            end

            obj.res = doLoadFiles(obj.hCfg);
            if isfield(obj.res, 'hClust')
                obj.res.hClust.hCfg = obj.hCfg;
            end
        end
    end

    % GETTERS/SETTERS
    methods
        % hCfg
        function set.hCfg(obj, hCfg)
            assert(isa(hCfg, 'jrclust.Config'), 'hCfg must be a Config object');
            obj.hCfg = hCfg;
        end
        
        % hClust
        function val = get.hClust(obj)
            if isempty(obj.res) || ~isfield(obj.res, 'hClust')
                val = [];
            else
                val = obj.res.hClust;
            end
        end

        % isError
        function val = get.isError(obj)
            val = obj.isError;
            if ~isempty(obj.hCfg) % pick up error in Config
                val = val || obj.hCfg.isError;
            end
            if ~isempty(obj.hDetect) % pick up error in hDetect
                val = val || obj.hDetect.isError;
            end
            if ~isempty(obj.hSort) % pick up error in hSort
                val = val || obj.hSort.isError;
            end
        end
    end
end
