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

        probe = loadProbe(pf);
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
            if ischar(userParams.(paramName)) && ~isempty(userParams.(paramName))
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