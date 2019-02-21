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
        clearMemory(obj);
        deprecateCmd(obj, oldCmd, iMsg, newCmd);
        error(obj, emsg);
        startParPool(obj);
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
