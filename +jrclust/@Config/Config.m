classdef Config < dynamicprops
    %CONFIG JRCLUST session configuration (P++)
    %% OBJECT-LEVEL PROPERTIES
    properties (Hidden, SetAccess=private, SetObservable)
        customParams;           % params not included in the default set
        isError;                % true if an error in configuration
        isV3Import;             % true if old-style params are referred to in configFile
        logEntries;             % containers.Map of log entries
        logFid;                 % file id of open log file
        paramSet;               % common and advanced parameter sets with default values and validation criteria
        tempParams;             % temporary parameters (probably a hack)
    end

    %% DEPENDENT OBJECT-LEVEL PROPERTIES
    properties (Dependent, Hidden, SetObservable)
        deprecatedParams;       % deprecated parameters (for excluding)
        fullParams;             % all parameters (merged common/advanced param sets)
        oldParamSet;            % old-style parameters with mapping to new params
    end

    %% CONFIG FILE
    properties (SetAccess=private, SetObservable)
        configFile;             % path to configuration file
    end

    %% NOT TO BE SET BY USERS
    properties (SetAccess=private)
        siteNeighbors;          % indices of neighbors for each site
    end

    %% RECORDING(S) (to ease the v3-v4 transition)
    properties (Dependent, Hidden, SetObservable)
        singleRaw;              % formerly vcFile
        multiRaw;               % formerly csFile_merge
    end

    %% COMPUTED PARAMS
    properties (SetObservable, Dependent)
        bytesPerSample;         % byte count for each discrete sample
        evtManualThreshSamp;    % evtManualThresh / bitScaling
        evtWindowRawSamp;       % interval around event to extract raw spike waveforms, in samples
        evtWindowSamp;          % interval around event to extract filtered spike waveforms, in samples
        nSites;                 % numel(siteMap)
        nSitesEvt;              % 2*nSiteDir + 1 - nSitesExcl
        refracIntSamp;          % spike refractory interval, in samples

        featuresFile;           % path to binary file containing spike features
        filtFile;               % path to binary file containing filtered traces
        rawFile;                % path to binary file containing raw traces
        resFile;                % path to MAT file containing results struct
        sessionName;            % name of prm file, without path or extensions
    end

    %% LIFECYCLE
    methods
        function obj = Config(filename)
            %CONFIG Construct an instance of this class
            skipValidation = 0;

            if nargin == 0
                userParams = struct();
                obj.configFile = '';
                skipValidation = 1;
            elseif isstruct(filename)
                userParams = filename;
                obj.configFile = '';
            elseif ischar(filename)
                filename_ = jrclust.utils.absPath(filename);
                if isempty(filename_)
                    error('File not found: ''%s''', filename);
                end

                userParams = jrclust.utils.mToStruct(filename_);
                obj.configFile = filename_;
            end

            obj.customParams = {};
            obj.isV3Import = 0;
            obj.isError = 0;

            obj.logEntries = containers.Map();

            % for setting temporary parameters
            obj.tempParams = containers.Map();

            obj.loadParams(userParams);
            if ~skipValidation
                obj.validateParams();
            end

            % define a default outputDir if not already set
            if ~isempty(obj.configFile) && isempty(obj.outputDir)
                obj.outputDir = fileparts(obj.configFile);
            elseif isempty(obj.outputDir)
                obj.outputDir = pwd();
            end
        end

        function obj = subsasgn(obj, prop, val)
            if strcmp(prop.type, '.')
                [flag, val, errMsg] = obj.validateProp(prop.subs, val);
                if flag
                    obj.setProp(prop.subs, val);
                else
                    error(errMsg);
                end
            end
        end
    end

    %% DOUBLE SECRET METHODS
    methods (Access = private, Hidden)
        loadParams(obj, filename);
        setCustomProp(obj, propname, val);
        setProp(obj, propname, val);
        validateParams(obj);
        [flag, val, errMsg] = validateProp(obj, propname, val);
        warning(obj, wmsg, varargin);
    end

    %% SECRET METHODS
    methods (Hidden)
        setConfigFile(obj, configFile, reloadParams);
    end

    %% USER METHODS
    methods
        function val = getOr(obj, fn, dv)
            %GETOR GET set value `obj.(fn)` OR default value `dv` if unset or empty
            if nargin < 3
                dv = [];
            end

            if ~isprop(obj, fn) || isempty(obj.(fn))
                val = dv;
            else
                val = obj.(fn);
            end
        end

        function rd = recDurationSec(obj, recID)
            %RECDURATIONSEC Get duration of recording file(s) in seconds
            if nargin < 2 || isempty(recID)
                hRecs = cellfun(@(fn) jrclust.detect.newRecording(fn, obj), obj.rawRecordings, 'UniformOutput', 0);
                rd = sum(cellfun(@(hR) hR.nSamples, hRecs))/obj.sampleRate;
            elseif recID < 1 || recID > numel(obj.rawRecordings)
                error('recording ID %d is invalid (there are %d recordings)', recID, numel(obj.rawRecordings));
            else
                hRec = jrclust.detect.newRecording(obj.rawRecordings{recID}, obj);
                rd = hRec.nSamples/obj.sampleRate;
            end
        end

        function resetTemporaryParams(obj, prmKeys)
            %RESETTEMPORARYPARAMS Reset temporary parameters
            if nargin < 2 || isempty(prmKeys)
                prmKeys = keys(obj.tempParams);
            elseif nargin == 2
                if ischar(prmKeys)
                    prmKeys = {prmKeys};
                end
                % only try to reset parameters we actually have
                prmKeys = intersect(prmKeys, keys(obj.tempParams));
            end

            for i = 1:numel(prmKeys)
                fn = prmKeys{i};
                obj.(fn) = obj.tempParams(fn);
                remove(obj.tempParams, fn);
            end
        end

        function setTemporaryParams(obj, varargin)
            %SETTEMPORARYPARAMS Set temporary parameters to reset later
            prmKeys = varargin(1:2:end);
            prmVals = varargin(2:2:end);

            if numel(prmKeys) ~= numel(prmVals)
                warning('number of property names not equal to values; skipping');
                return;
            end

            for i = 1:numel(prmKeys)
                prmKey = prmKeys{i};

                % already set a temporary value for this parameter, reset
                % it or we'll lose the original
                if isKey(obj.tempParams, prmKey)
                    obj.resetTemporaryParams(prmKey);
                end
                try
                    obj.tempParams(prmKey) = obj.(prmKey); % save old value for later
                    obj.(prmKey) = prmVals{i};
                catch ME
                    remove(obj.tempParams, prmKey);
                    warning('failed to set %s: %s', prmKey, ME.message);
                end
            end
        end
    end

    %% GETTERS/SETTERS
    methods
        % bytesPerSample
        function bp = get.bytesPerSample(obj)
            bp = jrclust.utils.typeBytes(obj.dataType);
        end

        % deprecatedParams
        function val = get.deprecatedParams(obj)
            if ~isstruct(obj.paramSet)
                val = [];
            else
                val = obj.paramSet.deprecated;
            end
        end

        % evtManualThreshSamp
        function mt = get.evtManualThreshSamp(obj)
            mt = obj.evtManualThresh / obj.bitScaling;
        end

        % evtWindowRawSamp
        function ew = get.evtWindowRawSamp(obj)
            if isprop(obj, 'evtWindowRaw') && isprop(obj, 'sampleRate')
                ew = round(obj.evtWindowRaw * obj.sampleRate / 1000);
            else
                ew = [];
            end
        end

        % evtWindowSamp
        function ew = get.evtWindowSamp(obj)
            if isprop(obj, 'evtWindow') && isprop(obj, 'sampleRate')
                ew = round(obj.evtWindow * obj.sampleRate / 1000);
            else
                ew = [];
            end
        end

        % featuresFile
        function val = get.featuresFile(obj)
            val = fullfile(obj.outputDir, [obj.sessionName '_features.jrc']);
        end

        % filtFile
        function val = get.filtFile(obj)
            val = fullfile(obj.outputDir, [obj.sessionName '_filt.jrc']);
        end

        % fullParams
        function val = get.fullParams(obj)
            val = jrclust.utils.mergeStructs(obj.paramSet.commonParameters, ...
                                             obj.paramSet.advancedParameters);
        end

        % isError
        function val = get.isError(obj)
            val = obj.isError;
            if isprop(obj, 'rawRecordings')
                switch class(obj.rawRecordings)
                    case 'cell'
                        val = val | all(cellfun(@isempty, obj.rawRecordings));

                    case 'char'
                        val = val | isempty(obj.rawRecordings);

                    otherwise
                        val = 1;
                end
            end
        end

        % multiRaw
        function set.multiRaw(obj, mr)
            if ~isprop(obj, 'rawRecordings')
                addprop(obj, 'rawRecordings');
            end
            if ~isprop(obj, 'configFile')
                addprop(obj, 'configFile');
                obj.configFile = '';
            end

            if ischar(mr) && ~any(mr == '*')
                obj.singleRaw = mr;
                return;
            elseif ischar(mr) % wildcard character
                if isprop(obj, 'configFile') && ~isempty(obj.configFile)
                    basedir = fileparts(obj.configFile);
                else
                    basedir = pwd();
                end

                mr_ = jrclust.utils.absPath(mr, basedir);
                if isempty(mr_)
                    error('Wildcard not recognized: %s', mr);
                end
            else
                % check is a cell array
                assert(iscell(mr), 'multiRaw must be a cell array');

                % get absolute paths
                if isprop(obj, 'configFile') && ~isempty(obj.configFile)
                    basedir = fileparts(obj.configFile);
                else
                    basedir = pwd();
                end

                mr_ = cellfun(@(fn) jrclust.utils.absPath(fn, basedir), mr, 'UniformOutput', 0);

                % handle glob expansions
                while any(cellfun(@iscell, mr_))
                    mr_ = [mr_{:}];
                end

                isFound = ~cellfun(@isempty, mr_);
                if ~all(isFound)
                    error('%d/%d files not found', sum(isFound), numel(isFound));
                end
            end

            % validation done, just set prop
            obj.setProp('rawRecordings', jrclust.utils.sortNat(unique(mr_)));
        end

        % nSites
        function ns = get.nSites(obj)
            if isprop(obj, 'siteMap')
                ns = numel(obj.siteMap);
            else
                ns = [];
            end
        end

        % nSitesEvt
        function ns = get.nSitesEvt(obj)
            if isprop(obj, 'nSiteDir') && isprop(obj, 'nSitesExcl')
                ns = 2*obj.nSiteDir - obj.nSitesExcl + 1;
            else
                ns = [];
            end
        end

        % oldParamSet
        function val = get.oldParamSet(obj)
            if ~isstruct(obj.paramSet)
                val = [];
            else
                val = obj.paramSet.old2new;
            end
        end

        % rawFile
        function val = get.rawFile(obj)
            val = fullfile(obj.outputDir, [obj.sessionName '_raw.jrc']);
        end

        % refracIntSamp
        function ri = get.refracIntSamp(obj)
            if isprop(obj, 'refracInt') && isprop(obj, 'sampleRate')
                ri = round(obj.refracInt * obj.sampleRate / 1000);
            else
                ri = [];
            end
        end

        % resFile
        function val = get.resFile(obj)
            val = fullfile(obj.outputDir, [obj.sessionName '_res.mat']);
        end

        % sessionName
        function sn = get.sessionName(obj)
            if isprop(obj, 'configFile') && ~isempty(obj.configFile)
                [~, sn, ~] = fileparts(obj.configFile);
            else
                sn = '';
            end
        end

        % singleRaw
        function set.singleRaw(obj, sr)
            if ~isprop(obj, 'rawRecordings')
                addprop(obj, 'rawRecordings');
            end
            if ~isprop(obj, 'configFile')
                addprop(obj, 'configFile');
                obj.configFile = '';
            end

            if iscell(sr) || (ischar(sr) && any(sr == '*'))
                obj.multiRaw = sr;
                return;
            end

            % check is a cell array
            assert(ischar(sr), 'singleRaw must be a string');

            % get absolute paths
            if isprop(obj, 'configFile') && ~isempty(obj.configFile)
                basedir = fileparts(obj.configFile);
            else
                basedir = pwd();
            end
            sr_ = jrclust.utils.absPath(sr, basedir);
            if isempty(sr_)
                error('Recording ''%s'' not found', sr);
            end

            % validation done, just set prop
            obj.setProp('rawRecordings', {sr_});
        end
    end
end
