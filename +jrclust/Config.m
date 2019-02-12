classdef Config < dynamicprops
    %CONFIG JRCLUST session configuration
    % replacement for P struct

    %% OBJECT-LEVEL PROPERTIES
    properties (Hidden, SetAccess=private, SetObservable)
        customParams;           % params not included in the default set
        isError;                % true if an error in configuration
        isV3Import;             % true if old-style params are referred to in configFile
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
        sessionName;            % name of prm file, without path or extensions
    end

    %% LIFECYCLE
    methods
        function obj = Config(filename)
            %CONFIG Construct an instance of this class
            if nargin == 0
                userParams = struct();
                obj.configFile = '';
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

            % for setting temporary parameters
            obj.tempParams = containers.Map();

            obj.loadParams(userParams);
            if ~isempty(obj.configFile) % prm file was specified, validate
                obj.validateParams();
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
        function error(obj, emsg, varargin)
            %ERROR Raise an error
            obj.isError = 1;
            if obj.batchMode
                error(emsg);
            else
                uiwait(errordlg(emsg, varargin{:}, 'modal'));
            end
        end

        function loadParams(obj, filename)
            %LOADPARAMS Load parameters from file
            if nargin < 2
                filename = obj.configFile;
            end

            if ischar(filename)
                filename_ = jrclust.utils.absPath(filename);
                userParams = jrclust.utils.mToStruct(filename_); % raises error if not a file
            elseif isstruct(filename)
                userParams = filename;
            else
                error('Class not recognized: %s', class(filename));
            end

            % read in default parameter set
            obj.paramSet = jrclust.utils.getDefaultParams(0);

            % set default parameters
            paramNames = fieldnames(obj.fullParams);
            for i = 1:numel(paramNames)
                paramName = paramNames{i};
                if strcmp(paramName, 'rawRecordings')
                    obj.setProp('rawRecordings', '');
                    continue;
                end
                [flag, val, errMsg] = obj.validateProp(paramName, obj.fullParams.(paramName).default_value);
                if flag
                    obj.setProp(paramName, val);
                else
                    warning(errMsg);
                end
            end

            % overwrite default parameters with user-specified params
            if isfield(userParams, 'template_file') && ~isempty(userParams.template_file)
                try
                    userTemplate = jrclust.utils.mToStruct(jrclust.utils.absPath(userParams.template_file));
                    fns = fieldnames(userParams);
                    for i = 1:numel(fns) % merge userParams (specific) into userTemplate (general)
                        userTemplate.(fns{i}) = userParams.(fns{i});
                    end
                    userParams = userTemplate;
                catch ME
                    obj.warning(sprintf('Could not set template file %s: %s ', userParams.template_file, ME.message), 'Missing template file');
                end
            end

            % set batchMode first because it is used in the loop
            if isfield(userParams, 'batchMode')
                [flag, val, errMsg] = obj.validateProp('batchMode', userParams.batchMode);
                if flag
                    obj.setProp('batchMode', val);
                else
                    warning(errMsg);
                end

                userParams = rmfield(userParams, 'batchMode');
            end

            % load probe from a probe file (legacy support)
            if isfield(userParams, 'probe_file') && ~isempty(userParams.probe_file)
                % first check local directory
                if isprop(obj, 'configFile') && ~isempty(obj.configFile)
                    basedir = fileparts(obj.configFile);
                else
                    basedir = fullfile(jrclust.utils.basedir(), 'probes');
                end

                pf = jrclust.utils.absPath(userParams.probe_file, basedir);
                if isempty(pf)
                    pf = jrclust.utils.absPath(userParams.probe_file, fullfile(jrclust.utils.basedir(), 'probes'));
                end
                if isempty(pf)
                    obj.error(sprintf('Could not find probe file ''%s''', userParams.probe_file));
                end

                probe = doLoadProbe(pf);
                probeFields = fieldnames(probe);

                for i = 1:numel(probeFields)
                    fn = probeFields{i};
                    obj.(fn) = probe.(fn);
                end

                userParams = rmfield(userParams, 'probe_file');
            end

            % set user-specified params
            uParamNames = fieldnames(userParams);
            for i = 1:numel(uParamNames)
                paramName = uParamNames{i};

                % ignore configFile/template_file
                if ismember(paramName, {'configFile', 'vcFile_prm', 'template_file'})
                    continue;
                elseif strcmpi(paramName, 'vcFile') && isempty(userParams.vcFile)
                    continue;
                elseif strcmpi(paramName, 'csFile_merge') && isempty(userParams.csFile_merge)
                    continue;
                end

                % empty values in the param file take on their defaults
                if strcmp(paramName, 'rawRecordings') % validate depending on type
                    if ischar(userParams.(paramName))
                        obj.singleRaw = userParams.(paramName);
                    elseif iscell(userParams.(paramName))
                        obj.multiRaw = userParams.(paramName);
                    else
                        obj.error('rawRecordings must be a char or cell array of char');
                    end

                    if isempty(obj.rawRecordings)
                        obj.error('rawRecordings cannot be empty');
                    end
                elseif ~isempty(userParams.(paramName))
                    [flag, val, errMsg] = obj.validateProp(paramName, userParams.(paramName));
                    if flag
                        obj.setProp(paramName, val);
                    else % TODO: warn users after a grace period
                        % warning(errMsg);
                    end
                end
            end
        end

        function setCustomProp(obj, propname, val)
            %SETCUSTOMPROP Set a property not included in the defaults
            if ismember(propname, obj.deprecatedParams.unsupported)
                return;
            end

            % ignore property if it's Dependent
            propData = ?jrclust.Config;
            propNames = {propData.PropertyList.Name};
            dependentProps = propNames([propData.PropertyList.Dependent]);
            if ismember(propname, dependentProps)
                return;
            end

            if ~isprop(obj, propname)
                addprop(obj, propname);
                if ~ismember(propname, obj.customParams)
                    obj.customParams{end+1} = propname;
                end
            end

            obj.(propname) = val;
        end

        function setProp(obj, propname, val)
            %SETPROP Set a property
            if isfield(obj.oldParamSet, propname)
                propname = obj.oldParamSet.(propname);
                obj.isV3Import = 1;
            end

            if isfield(obj.fullParams, propname)
                if ~isprop(obj, propname)
                    obj.addprop(propname);
                end

                obj.(propname) = val;
            elseif ismember(propname, {'singleRaw', 'multiRaw'}) % separate validation for these
                obj.(propname) = val;
            else
                obj.setCustomProp(propname, val);
            end
        end

        function validateParams(obj)
            %VALIDATEPARAMS Validate parameters and compute others
            if obj.nSites == 0
                obj.error('No siteMap specified', 'Bad probe configuration');
            end

            if size(obj.siteLoc, 1) ~= obj.nSites
                obj.error('Malformed probe geometry', 'Bad probe configuration');
                return;
            end

            if numel(obj.shankMap) ~= obj.nSites
                obj.error('Malformed shank indexing', 'Bad probe configuration');
                return;
            end

            if max(obj.siteMap) > obj.nChans
                obj.error('siteMap refers to channels larger than indexed by nChans', 'Bad probe configuration');
                return;
            end

            % nSiteDir and/or nSitesExcl may not have been specified
            if isempty(obj.nSiteDir) || isempty(obj.nSitesExcl)
                siteDists = pdist2(obj.siteLoc, obj.siteLoc);

                % max over all sites of number of neighbors in detect radius
                nNeighDetect = max(sum(siteDists <= obj.evtDetectRad)); % 11/7/17 JJJ: med to max

                if isempty(obj.nSitesExcl)
                    % max over all sites of number of neighbors in extract radius
                    nNeighExtract = max(sum(siteDists <= obj.evtGroupRad)); % 11/7/17 JJJ: med to max
                    nsd = (nNeighExtract - 1)/2;
                    obj.nSitesExcl = nNeighExtract - nNeighDetect;
                else
                    nNeighExtract = nNeighDetect + obj.nSitesExcl;
                    nsd = (nNeighExtract - 1)/2;
                end

                if isempty(obj.nSiteDir)
                    obj.nSiteDir = nsd;
                end
            end

            if obj.nSitesEvt <= 0
                obj.error('nSitesExcl is too large or nSiteDir is too small', 'Bad configuration');
            end

            % ignoreSites/ignoreChans
            obj.ignoreChans = obj.ignoreChans(ismember(obj.ignoreChans, obj.siteMap));
            obj.ignoreSites = intersect(obj.ignoreSites, 1:numel(obj.siteMap));
            obj.ignoreSites = union(obj.ignoreSites, find(ismember(obj.siteMap, obj.ignoreChans)));

            obj.siteNeighbors = findSiteNeighbors(obj.siteLoc, 2*obj.nSiteDir + 1, obj.ignoreSites, obj.shankMap);

            % boost that gain
            obj.bitScaling = obj.bitScaling/obj.gainBoost;
        end

        function [flag, val, errMsg] = validateProp(obj, propname, val)
            %VALIDATEPROP Ensure a property is valid
            if isfield(obj.oldParamSet, propname) % map the old param name to the new one
                propname = obj.oldParamSet.(propname);
                obj.isV3Import = 1;
            end

            flag = 1;
            errMsg = '';

            % found in default params, do validation
            if isfield(obj.fullParams, propname)
                validData = obj.fullParams.(propname).validation;
                classes = validData.classes;
                attributes = validData.attributes;

                if isempty(val) || isempty(attributes)
                    if ~any(cellfun(@(c) isa(val, c), classes))
                        flag = 0;
                    end

                    return;
                end

                try
                    validateattributes(val, classes, attributes);

                    % transform val in some way
                    if isfield(validData, 'postapply')
                        hFun = eval(validData.postapply);
                        val = hFun(val);
                    end

                    % check additional constraints
                    if isfield(validData, 'postassert')
                        hFun = eval(validData.postassert);
                        assert(hFun(val));
                    end

                    if isfield(validData, 'values')
                        assert(all(ismember(val, validData.values)));
                    end
                catch ME
                    errMsg = sprintf('Could not set %s: %s', propname, ME.message);
                    flag = 0;
                end
            end
        end

        function warning(obj, wmsg, varargin)
            %WARNING Raise a warning
            if obj.batchMode
                warning(wmsg);
            else
                uiwait(warndlg(wmsg, varargin{:}, 'modal'));
            end
        end
    end

    %% SECRET METHODS
    methods (Hidden)
        function setConfigFile(obj, configFile, reloadParams)
            %SETCONFIGFILE Don't use this. Seriously.
             if nargin < 2
                 return;
             end
             if nargin < 3
                 reloadParams = 1;
             end

             configFile_ = jrclust.utils.absPath(configFile);
             if isempty(configFile_)
                 error('Could not find %s', configFile);
             end

             obj.configFile = configFile_;

             if reloadParams
                 obj.loadParams(obj.configFile);
             end
        end
    end

    %% USER METHODS
    methods
        function edit(obj)
            %EDIT Edit the config file
            edit(obj.configFile);
        end

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

        function success = save(obj, filename, exportAdv, diffsOnly)
            %SAVE Write parameters to a file
            success = 0;

            if nargin < 2
                filename = obj.configFile;
            end
            if nargin < 3
                exportAdv = 0;
            end
            if nargin < 4
                diffsOnly = 0;
            end

            if isempty(filename) % passed an empty string or no config file
                filename = 'stdout';
            end

            if ~strcmpi(filename, 'stdout')
                filename_ = jrclust.utils.absPath(filename);
                if isempty(filename_)
                    error('Could not find ''%s''', filename);
                elseif exist(filename, 'dir')
                    error('''%s'' is a directory', filename);
                end

                filename = filename_;
            end

            if strcmpi(filename, 'stdout')
                fid = 1;
            else
                if isempty(obj.configFile) % bind configFile to this new path
                    obj.configFile = filename;
                end

                % file already exists, back it up!
                if exist(filename, 'file')
                    [~, ~, ext] = fileparts(filename);
                    backupFile = jrclust.utils.subsExt(filename, [ext, '.bak']);
                    try
                        copyfile(filename, backupFile);
                    catch ME % cowardly back out
                        warning('Could not back up old file: %s', ME.message);
                        return;
                    end
                end

                [fid, errmsg] = fopen(filename, 'w');
                if fid == -1
                    warning('Could not open config file for writing: %s', errmsg);
                    return;
                end
            end

            paramsToExport = obj.paramSet.commonParameters;
            if exportAdv
                paramsToExport = jrclust.utils.mergeStructs(paramsToExport, obj.paramSet.advancedParameters);
            end

            % replace fields in paramsToExport with values in this object
            paramNames = fieldnames(paramsToExport);
            for i = 1:numel(paramNames)
                pn = paramNames{i};

                if jrclust.utils.isEqual(paramsToExport.(pn).default_value, obj.(pn))
                    if diffsOnly % don't export fields which have default values
                        paramsToExport = rmfield(paramsToExport, pn);
                    end
                else
                    paramsToExport.(pn).default_value = obj.(pn);
                end
            end

            % write the file
            paramNames = fieldnames(paramsToExport);
            sections = {'usage', 'execution', 'probe', 'recording file', ...
                        'preprocessing', 'spike detection', 'feature extraction', ...
                        'clustering', 'curation', 'display', 'trial', ...
                        'validation', 'preview', 'traces', 'lfp', 'aux channel'};
            [~, new2old] = jrclust.utils.getOldParamMapping();

            % write header
            progInfo = jrclust.utils.info();
            fprintf(fid, '%% %s parameters ', progInfo.program);
            if ~exportAdv
                fprintf(fid, '(common parameters only) ');
            end
            if diffsOnly
                fprintf(fid, '(default parameters not exported)');
            end
            fprintf(fid, '\n%% For a description of these parameters, see %s\n\n', [progInfo.docsSite, 'parameters/index.html']);

            % write sections
            for i = 1:numel(sections)
                section = sections{i};
                % no params have this section as primary, skip it
                if ~any(cellfun(@(pn) strcmp(section, paramsToExport.(pn).section{1}), paramNames))
                    continue;
                end

                fprintf(fid, '%% %s PARAMETERS\n', upper(section));

                for j = 1:numel(paramNames)
                    pn = paramNames{j};
                    pdata = paramsToExport.(pn);
                    if ~strcmpi(pdata.section{1}, section)
                        continue;
                    end

                    fprintf(fid, '%s = %s; %% ', pn, jrclust.utils.field2str(obj.(pn)));
                    if isfield(new2old, pn) % write old parameter name
                        fprintf(fid, '(formerly %s) ', new2old.(pn));
                    end
                    fprintf(fid, '%s', strrep(pdata.description, 'μ', char(956))); % \mu
                    if isempty(pdata.comment)
                        fprintf(fid, '\n');
                    else
                        fprintf(fid, ' (%s)\n', strrep(pdata.comment, 'μ', char(956))); % \mu
                    end
                end

                fprintf(fid, '\n');
            end

            % write out custom parameters
            if ~isempty(obj.customParams)
                fprintf(fid, '%% USER-DEFINED PARAMETERS\n');
                for j = 1:numel(obj.customParams)
                    pn = obj.customParams{j};

                    fprintf(fid, '%s = %s;', pn, jrclust.utils.field2str(obj.(pn)));
                    if isfield(new2old, pn) % write old parameter name
                        fprintf(fid, ' %% (formerly %s)\n', new2old.(pn));
                    else
                        fprintf(fid, '\n');
                    end
                end
            end

            if fid > 1
                fclose(fid);
            end

            success = 1;
        end

        function rd = recDurationSec(obj, recID)
            %RECDURATIONSECS Get duration of recording file(s) in seconds
            if nargin < 2 || isempty(recID)
                hRecs = cellfun(@(fn) jrclust.models.recording.Recording(fn, obj), obj.rawRecordings, 'UniformOutput', 0);
                rd = sum(cellfun(@(hR) hR.nSamples, hRecs))/obj.sampleRate;
            elseif recID < 1 || recID > numel(obj.rawRecordings)
                error('recording ID %d is invalid (there are %d recordings)', recID, numel(obj.rawRecordings));
            else
                hRec = jrclust.models.recording.Recording(obj.rawRecordings{recID}, obj);
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

        % fullParams
        function val = get.fullParams(obj)
            val = jrclust.utils.mergeStructs(obj.paramSet.commonParameters, ...
                                             obj.paramSet.advancedParameters);
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

        % refracIntSamp
        function ri = get.refracIntSamp(obj)
            if isprop(obj, 'refracInt') && isprop(obj, 'sampleRate')
                ri = round(obj.refracInt * obj.sampleRate / 1000);
            else
                ri = [];
            end
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
                error('''%s'' not found', sr);
            end

            % validation done, just set prop
            obj.setProp('rawRecordings', {sr_});
        end
    end
end
