classdef JRC < handle & dynamicprops
    %JRC

    properties (SetObservable, SetAccess=private, Hidden)
        errMsg;
        isCompleted;
        isError;
        isDetection;
        isSorting;
        isCuration;
    end

    properties (SetObservable, SetAccess=private)
        args;
        cmd;
        hCfg;
        hDet;
        hSort;
        hCurate;
    end

    % most relevant data from detection step
    properties (Dependent)
        spikeTimes;
        spikeSites;
    end

    % contains all computed data from detection step
    properties (Hidden, SetAccess=private)
        dRes;
    end

    % LIFECYCLE
    methods
        function obj = JRC(varargin)
            %JRC Construct an instance of this class
            obj.args = varargin;
            obj.isCompleted = false;
            obj.isError = false;
            obj.isDetection = false;
            obj.isSorting = false;
            obj.isCuration = false;

            if ~jrclust.utils.sysCheck()
                obj.errMsg = 'system requirements not met';
                obj.isError = true;
            elseif nargin > 0
                % handle arguments (legacy mode)
                obj.processArgs();
            end
        end
    end

    methods (Access=protected)
        function processArgs(obj)
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
                    imsg = 'Use `hJRC.hCfg = myconfig;` instead';
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
            legalCmds = {'detect', 'full', 'manual', 'sort', 'spikesort'};
            if ~any(strcmpi(obj.cmd, legalCmds))
                obj.errMsg = sprintf('Command `%s` not recognized', obj.cmd);
                errordlg(obj.errMsg, 'Unrecognized command');
                obj.isError = true;
                return;
            end

            % determine which commands in the pipeline to run
            detectCmds = {'detect', 'sort', 'spikesort', 'full'};
            sortCmds   = {'sort', 'spikesort', 'full'};
            curateCmds = {'full', 'manual'};

            if any(strcmp(obj.cmd, detectCmds))
                obj.isDetection = true;
            end
            if any(strcmp(obj.cmd, sortCmds))
                obj.isSorting = true;
            end
            if any(strcmp(obj.cmd, curateCmds))
                obj.isCuration = true;
            end

            % load from saved
            % TODO (much much later)

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

    % USER METHODS
    methods
        function ip = inProgress(obj)
            %INPROGRESS 
            ip = ~(obj.isCompleted || obj.isError);
        end

        function it = invocationType(obj)
            if obj.isError
                it = 'error';
            elseif obj.isSorting
                it = 'sorting';
            elseif obj.isCuration
                it = 'manual';
            else % ~(obj.isCuration || obj.isSorting)
                it = 'info';
            end
        end

        function run(obj)
            if obj.isError
                error(obj.errMsg);
            elseif obj.isCompleted
                warning('command %s completed successfully; to rerun, use rerun()');
                return;
            end

            if obj.isDetection
                obj.hDet = jrclust.controllers.DetectController(obj.hCfg);
                obj.dRes = obj.hDet.detect();
                if obj.hCfg.fVerbose
                    fprintf('Detection completed in %0.2f seconds', obj.dRes.runtime);
                end

                if ~(obj.isSorting || obj.isCuration) % save to pick up later
                    obj.hDet.saveFiles();
                end
            end

            obj.isCompleted = true;
        end

        function rerun(obj)
            if obj.isError
                error(obj.errMsg);
            else
                obj.isCompleted = false;
                obj.run();
            end
        end
    end

    % GETTERS/SETTERS
    methods
        % args
        function args = get.args(obj)
            args = obj.args;
        end
        function set.args(obj, args)
            obj.args = args;
        end

        % hCfg
        function hCfg = get.hCfg(obj)
            hCfg = obj.hCfg;
        end
        function set.hCfg(obj, hCfg)
            obj.hCfg = hCfg;
        end

        % isError
        function ie = get.isError(obj)
            ie = obj.isError || obj.hCfg.isError;
        end

        % spikeTimes
        function st = get.spikeTimes(obj)
            if isempty(obj.dRes) || ~isfield(obj.dRes, 'spikeTimes')
                st = [];
            else
                st = obj.dRes.spikeTimes;
            end
        end

        % spikeSites
        function ss = get.spikeSites(obj)
            if isempty(obj.dRes) || ~isfield(obj.dRes, 'spikeSites')
                ss = [];
            else
                ss = obj.dRes.spikeSites;
            end
        end
    end
end

