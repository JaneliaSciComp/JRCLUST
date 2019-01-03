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
        hDet;           %
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
        dRes;           % detect step (must contain at a minimum spikeTimes and spikeSites)
        sRes;           % sort step
        cRes;           % curate step
    end

    %% LIFECYCLE
    methods
        function obj = JRC(varargin)
            %JRC Construct an instance of this class
            obj.args = varargin;
            obj.isCompleted = false;
            obj.isError = false;
            obj.isDetect = false;
            obj.isSort = false;
            obj.isCurate = false;

            if ~jrclust.utils.sysCheck()
                obj.errMsg = 'system requirements not met';
                obj.isError = true;
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

            switch obj.cmd
                % deprecated commands; may be removed in a future release
                case {'compile-ksort', 'edit', 'git-pull', 'issue', 'import-kilosort-sort', ...
                      'import-ksort-sort', 'kilosort', 'kilosort-verify', 'ksort', 'ksort-verify' ...
                      'which', 'wiki', 'wiki-download'}
                    jrclust.utils.depWarn(obj.cmd);
                    obj.isCompleted = true;
                    return;

                case {'doc', 'doc-edit'}
                    imsg = 'Please visit the wiki at https://github.com/JaneliaSciComp/JRCLUST/wiki';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'download'
                    imsg = 'You can find sample.bin and sample.meta at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'gui'
                    imsg = 'GUI is not implemented yet, but eventually you can just use `jrc`';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'install'
                    imsg = 'You might be looking for `compile` instead';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case {'set', 'setprm', 'set-prm'}
                    imsg = 'Create a new JRC handle instead';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'update'
                    imsg = 'Please check the repository at https://github.com/JaneliaSciComp/JRCLUST for updates';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.isCompleted = true;
                    return;

                % deprecated synonyms, warn but proceed
                case 'spikedetect'
                    imsg = 'Please use ''detect'' in the future';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.cmd = 'detect';

                case {'cluster', 'clust', 'sort-verify', 'sort-validate', 'sort-manual'}
                    imsg = 'Please use ''sort'' in the future';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.cmd = 'sort';

                case {'detectsort', 'detect-sort', 'spikesort-verify', ...
                      'spikesort-validate', 'spikesort-manual', 'detectsort-manual'}
                    imsg = 'Please use ''spikesort'' in the future';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.cmd = 'spikesort';

                case 'all'
                    imsg = 'Please use ''full'' in the future';
                    jrclust.utils.depWarn(obj.cmd, imsg);
                    obj.cmd = 'full';

                % info commands
                case 'about'
                    md = jrclust.utils.info();
                    verstr = sprintf('%s v%s', md.program, jrclust.utils.version());
                    abstr = jrclust.utils.about();
                    msgbox(abstr, verstr);
                    obj.isCompleted = true;
                    return;

                case 'help'
                    disp(jrclust.utils.help());
                    obj.isCompleted = true;
                    return;

                case 'version'
                    md = jrclust.utils.info();
                    fprintf('%s v%s\n', md.program, jrclust.utils.version());
                    obj.isCompleted = true;
                    return;
            end

            % command sentinel
            legalCmds = {'detect', 'sort', 'manual', 'full'};
            if ~any(strcmpi(obj.cmd, legalCmds))
                obj.errMsg = sprintf('Command `%s` not recognized', obj.cmd);
                errordlg(obj.errMsg, 'Unrecognized command');
                obj.isError = true;
                return;
            end

            % determine which commands in the pipeline to run
            detectCmds = {'detect', 'full'};
            sortCmds   = {'sort', 'spikesort', 'full'};
            curateCmds = {'manual', 'full'};

            if any(strcmp(obj.cmd, curateCmds))
                obj.isCurate = true;
            end

            if any(strcmp(obj.cmd, sortCmds))
                obj.isSort = true;
            end

            if any(strcmp(obj.cmd, detectCmds))
                obj.isDetect = true;
            end

            % commands from here on out require a parameter file
            if nargs < 1
                obj.errMsg = sprintf('Command `%s` requires a parameter file', obj.cmd);
                errordlg(obj.errMsg, 'Missing parameter file');
                obj.isError = true;
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

        function run(obj)
            %RUN Run commands
            if obj.isError
                error(obj.errMsg);
            elseif obj.isCompleted
                warning('command ''%s'' completed successfully; to rerun, use rerun()', obj.cmd);
                return;
            end

            % try to load sort and detect results
            if obj.isCurate && ~obj.isSort
                [dRes_, sRes_] = obj.loadFiles();
                if isempty(dRes_)
                    obj.isDetect = true;
                    obj.isSort = true;
                elseif isempty(sRes_)
                    obj.dRes = dRes_;
                    obj.isSort = true;
                else
                    obj.dRes = dRes_;
                    obj.sRes = sRes_;
                end
            end

            % try to load detect results
            if obj.isSort && ~obj.isDetect
                dRes_ = obj.loadFiles();
                if isempty(dRes_)
                    obj.isDetect = true;
                else
                    obj.dRes = dRes_;
                end
            end

            % notify user that we're using previously-computed results
            if obj.hCfg.verbose
                if ~isempty(obj.dRes) && isfield(obj.dRes, 'completedAt')
                    fprintf('Using spikes detected on %s\n', datestr(obj.dRes.completedAt));
                elseif ~isempty(obj.dRes)
                    fprintf('Using previously-detected spikes\n');
                end
                if ~isempty(obj.sRes) && isfield(obj.sRes, 'completedAt')
                    fprintf('Using clustering computed on %s\n', datestr(obj.sRes.completedAt));
                elseif ~isempty(obj.sRes)
                    fprintf('Using previously-clustered spikes\n');
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

            % DETECT SPIKES
            gpuDetect = obj.hCfg.useGPU; % save this in case useGPU is disabled during detection step
            if obj.isDetect
                obj.hDet = jrclust.controllers.detect.DetectController(obj.hCfg);
                obj.dRes = obj.hDet.detect();

                if obj.hDet.isError
                    error(obj.hDet.errMsg);
                elseif obj.hCfg.verbose
                    fprintf('Detection completed in %0.2f seconds\n', obj.dRes.runtime);
                end
            end

            % CLUSTER SPIKES
            gpuSort = obj.hCfg.useGPU | gpuDetect;
            if obj.isSort
                obj.hCfg.useGPU = gpuSort;

                obj.hSort = jrclust.controllers.sort.SortController(obj.hCfg);
                obj.sRes = obj.hSort.sort(obj.dRes);

                if obj.hSort.isError
                    error(obj.hSort.errMsg);
                elseif obj.hCfg.verbose
                    fprintf('Sorting completed in %0.2f seconds\n', obj.sRes.runtime);
                end
            end

            % CURATE SPIKES
            gpuCurate = obj.hCfg.useGPU | gpuSort;
            if obj.isCurate
                obj.hCfg.useGPU = gpuCurate;

                obj.hCurate = jrclust.controllers.curate.CurateController(obj.hClust);
                obj.hCurate.beginSession();
            end

            % save our results for later
            obj.saveFiles();
            obj.isCompleted = true;
        end

        function rerun(obj)
            %RERUN Rerun commands
            if obj.isError
                error(obj.errMsg);
            else
                obj.isCompleted = false;
                obj.run();
            end
        end

        function saveFiles(obj)
            %SAVEFILES Save results structs and binary files to disk
            if obj.isError
                error(obj.errMsg);
            end

            if isempty(obj.hCfg.outputDir)
                obj.hCfg.outputDir = fileparts(obj.hCfg.configFile);
            end

            doSaveFiles(obj);
        end

        function [dRes, sRes, cRes] = loadFiles(obj)
            %LOADFILES Load results structs
            if obj.isError
                error(obj.errMsg);
            end

            if isempty(obj.hCfg.outputDir)
                obj.hCfg.outputDir = fileparts(obj.hCfg.configFile);
            end

            [dRes, sRes, cRes] = deal([]);

            if nargout == 1
                dRes = doLoadFiles(obj);
            elseif nargout == 2
                [dRes, sRes] = doLoadFiles(obj);
            elseif nargout == 3
                [dRes, sRes, cRes] = doLoadFiles(obj);
            end
        end
    end

    % GETTERS/SETTERS
    methods
        % dRes
        function set.dRes(obj, dr)
            failMsg = 'dRes must contain at a minimum: spikeTimes, spikeSites';
            assert(isstruct(dr) && isfield(dr, 'spikeTimes') && isfield(dr, 'spikeSites'), failMsg);
            obj.dRes = dr;
        end

        % hClust
        function hc = get.hClust(obj)
            if isempty(obj.sRes)
                hc = [];
            else
                hc = obj.sRes.hClust;
            end
        end

        % isError
        function ie = get.isError(obj)
            ie = obj.isError;
            if ~isempty(obj.hCfg) % pick up error in Config
                ie = ie || obj.hCfg.isError;
            end
            if ~isempty(obj.hDet) % pick up error in hDet
                ie = ie || obj.hDet.isError;
            end
            if ~isempty(obj.hSort) % pick up error in hSort
                ie = ie || obj.hSort.isError;
            end
        end

        % spikeTimes
        function st = get.spikeTimes(obj)
            if isempty(obj.dRes)
                st = [];
            else
                st = obj.dRes.spikeTimes;
            end
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            if isempty(obj.dRes)
                ss = [];
            else
                ss = obj.dRes.spikeSites;
            end
        end

        % sRes
        function set.sRes(obj, sr)
            failMsg = 'sRes must contain at a minimum: hClust';
            assert(isstruct(sr) && isfield(sr, 'hClust') && isa(sr.hClust, 'jrclust.models.clustering.DensityPeakClustering'), failMsg);
            obj.sRes = sr;
        end
    end
end

