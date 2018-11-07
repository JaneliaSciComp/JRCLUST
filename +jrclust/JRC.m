classdef JRC < handle & dynamicprops
    %JRC
    
    properties
        config;
        args;
        doRun;
    end
    
    % lifecycle
    methods
        function obj = JRC(varargin)
            %JRC Construct an instance of this class
            jrclust.utils.sysCheck();
            obj.args = varargin;
            obj.doRun = false;

            % handle arguments (legacy mode)
            obj.processArgs();
        end

        function processArgs(obj)
            nArgs = numel(obj.args);

            if nArgs == 0
                return;
            end

            cmd = obj.args{1};
            disp(['your command was: ', cmd]);

            switch lower(cmd)
                % deprecated commands; may be removed in a future release
                case {'git-pull', 'issue', 'wiki', 'wiki-download', 'which', ...
                      'compile-ksort', 'kilosort', 'ksort', 'import-kilosort-sort', ...
                      'import-ksort-sort', 'kilosort-verify', 'ksort-verify'}
                    jrclust.utils.depWarn(cmd);
                    return;

                case 'download'
                    imsg = 'You can find sample.bin and sample.meta at https://drive.google.com/drive/folders/1-UTasZWB0TwFFFV49jSrpRPHmtve34O0?usp=sharing';
                    jrclust.utils.depWarn(cmd, imsg);
                    return;

                case {'doc', 'doc-edit'}
                    imsg = 'Please visit the wiki at https://github.com/JaneliaSciComp/JRCLUST/wiki';
                    jrclust.utils.depWarn(cmd, imsg);
                    return;

                case 'install'
                    imsg = 'You might be looking for `compile` instead';
                    jrclust.utils.depWarn(cmd, imsg);
                    return;

                case 'update'
                    imsg = 'Please check the repository at https://github.com/JaneliaSciComp/JRCLUST for updates';
                    jrclust.utils.depWarn(cmd, imsg);
                    return;

                case {'set', 'setprm', 'set-prm'}
                    imsg = 'Use `hJRC.config = myconfig` instead';
                    jrclust.utils.depWarn(cmd, imsg);
                    return;
            end
        end
    end
end

