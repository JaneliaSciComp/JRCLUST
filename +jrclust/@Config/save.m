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