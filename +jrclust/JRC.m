classdef JRC < handle & dynamicprops
    %JRC

    properties (SetObservable, SetAccess=private, Hidden)
        isCompleted;
        isError;
        isManualCuration;
        isSorting;
    end

    properties (SetObservable)
        args;
        hCfg;
    end

    % lifecycle
    methods
        function obj = JRC(varargin)
            %JRC Construct an instance of this class
            obj.args = varargin;
            obj.isCompleted = false;
            obj.isError = false;
            obj.isManualCuration = false;
            obj.isSorting = false;

            if ~jrclust.utils.sysCheck()
                obj.isError = true;
            elseif nargin > 0
                % handle arguments (legacy mode)
                obj.processArgs();
            end
        end

        function processArgs(obj)
            nArgs = numel(obj.args);

            if nArgs == 0
                return;
            end

            cmd = obj.args{1};
            switch lower(cmd)
                % deprecated commands; may be removed in a future release
                case {'compile-ksort', 'edit', 'git-pull', 'issue', 'import-kilosort-sort', ...
                      'import-ksort-sort', 'kilosort', 'kilosort-verify', 'ksort', 'ksort-verify' ...
                      'which', 'wiki', 'wiki-download'}
                    jrclust.utils.depWarn(cmd);
                    obj.isCompleted = true;
                    return;
                    
                case {'doc', 'doc-edit'}
                    imsg = 'Please visit the wiki at https://github.com/JaneliaSciComp/JRCLUST/wiki';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'download'
                    imsg = 'You can find sample.bin and sample.meta at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;
                    
                case 'gui'
                    imsg = 'GUI is not implemented yet, but eventually you can just use `jrc`';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'install'
                    imsg = 'You might be looking for `compile` instead';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;
                    
                case {'set', 'setprm', 'set-prm'}
                    imsg = 'Use `hJRC.hCfg = myconfig;` instead';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;

                case 'update'
                    imsg = 'Please check the repository at https://github.com/JaneliaSciComp/JRCLUST for updates';
                    jrclust.utils.depWarn(cmd, imsg);
                    obj.isCompleted = true;
                    return;

                % info commands
                case 'about'
                    disp(jrclust.utils.about());
                    obj.isCompleted = true;
                    return;

                case 'help'
                    disp(jrclust.utils.help());
                    obj.isCompleted = true;
                    return;

                case 'version'
                    md = jrclust.utils.info();
                    fprintf('%s %s\n', md.program, jrclust.utils.version());
                    obj.isCompleted = true;
                    return;
                    
            end

            % command sentinel
            legalCmds = {'sort'};
            if ~any(strcmpi(cmd, legalCmds))
                emsg = sprintf('Command `%s` not recognized', cmd);
                errordlg(emsg, 'Unrecognized command');
                obj.isError = true;
                return;
            end
            
            % load from saved
            % TODO (much much later)

            % commands from here on out require a parameter file
            if nArgs < 2
                emsg = sprintf('Command `%s` requires a parameter file', cmd);
                errordlg(emsg, 'Missing parameter file');
                obj.isError = true;
                return;
            end

            % load parameter file
            configFile = obj.args{2};
            obj.hCfg = jrclust.config.Config(configFile);
        end
    end
    
    methods % introspection
        function ip = inProgress(obj)
            ip = ~(obj.isCompleted || obj.isError);
        end

        function it = invocationType(obj)
            if obj.isError
                it = 'error';
            elseif obj.isSorting
                it = 'sorting';
            elseif obj.isManualCuration
                it = 'manual';
            else % ~(obj.isManualCuration || obj.isSorting)
                it = 'info';
            end
        end
    end

    methods % getters/setters
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
    end
end

