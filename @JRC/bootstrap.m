function bootstrap(obj, varargin)
    %DOBOOTSTRAP Bootstrap a JRCLUST session
    %   metafile: optional string; path (or glob) to meta file(s)
    if nargin > 1
        metafile_ = jrclust.utils.absPath(varargin{1});
        if isempty(metafile_) % warn?
            metafile = '';
            workingdir = pwd();
        elseif ischar(metafile_)
            workingdir = fileparts(metafile_);
            metafile = {metafile_};
        else % cell
            metafile = metafile_;
            workingdir = fileparts(metafile_{1});
        end

        if ~isempty(metafile)
            binfile = cellfun(@(f) jrclust.utils.subsExt(f, '.bin'), metafile, 'UniformOutput', 0);
        end
    else
        metafile = '';
        workingdir = pwd();
    end

    % first check for a .meta file
    if isempty(metafile)
        dlgAns = questdlg('Do you have a .meta file?', 'Bootstrap', 'No');

        switch dlgAns
            case 'Yes' % select .meta file
                [metafile, workingdir] = jrclust.utils.selectFile({'*.meta', 'SpikeGLX meta files (*.meta)'; '*.*', 'All Files (*.*)'}, 'Select one or more .meta files', workingdir, 1);
                if all(cellfun(@isempty, metafile))
                    return;
                end

                binfile = cellfun(@(f) jrclust.utils.subsExt(f, '.bin'), metafile, 'UniformOutput', 0);

            case 'No' % select recording file
                [binfile, workingdir] = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select one or more raw recordings', workingdir, 1);
                if all(cellfun(@isempty, binfile))
                    return;
                end

            case {'Cancel', ''}
                return;
        end
    end

    % check for missing binary files
    if any(cellfun(@(f) isempty(jrclust.utils.absPath(f)), binfile))
        binfile = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select one or more raw recordings', workingdir, 1);
        if cellfun(@isempty, binfile)
            return;
        end
    end

    % load metafile
    if ~isempty(metafile)
        SMeta_ = jrclust.utils.loadMetadata(metafile{1});
        cfgData = struct('sampleRate', SMeta_.sampleRate, ...
                         'nChans', SMeta_.nChans, ...
                         'bitScaling', SMeta_.bitScaling, ...
                         'headerOffset', 0, ...
                         'dataType', SMeta_.dataType, ...
                         'probe_file', fullfile(jrclust.utils.basedir(), 'probes', sprintf('%s.prb', SMeta_.probe)));

        cfgData.rawRecordings = binfile;
        cfgData.outputDir = workingdir;
    else
        cfgData.rawRecordings = binfile;
        cfgData.outputDir = workingdir;
    end

    dlgAns = questdlg('Would you like to specify a probe file?', 'Bootstrap', 'No');
    switch dlgAns
        case 'Yes' % select .prb file
            probedir = workingdir;
            if isempty(dir(fullfile(workingdir, '*.prb')))
                probedir = fullfile(jrclust.utils.basedir(), 'probes');
            end
            [probefile, probedir] = jrclust.utils.selectFile({'*.prb', 'Probe files (*.prb)'; '*.*', 'All Files (*.*)'}, 'Select a probe file', probedir, 0);
            cfgData.probe_file = fullfile(probedir, probefile);

        case {'Cancel', ''}
            return;
    end

    % construct the Config object from specified data
    hCfg_ = jrclust.Config(cfgData);

    while 1
        % confirm with the user
        [~, sessionName, ~] = fileparts(hCfg_.rawRecordings{1});
        configFile = fullfile(hCfg_.outputDir, [sessionName, '.prm']);

        dlgFieldNames = {'Config filename', ...
                         'Raw recording file(s)', ...
                         'Sampling rate (Hz)', ...
                         'Number of channels in file', ...
                         sprintf('%sV/bit', char(956)), ...
                         'Header offset (bytes)', ...
                         'Data Type (int16, uint16, single, double)'};
        dlgFieldVals = {configFile, ...
                        strjoin(hCfg_.rawRecordings, ','), ...
                        num2str(hCfg_.sampleRate), ...
                        num2str(hCfg_.nChans), ...
                        num2str(hCfg_.bitScaling), ...
                        num2str(hCfg_.headerOffset), ...
                        hCfg_.dataType};
        dlgAns = inputdlg(dlgFieldNames, 'Does this look correct?', 1, dlgFieldVals, struct('Resize', 'on', 'Interpreter', 'tex'));
        if isempty(dlgAns)
            return;
        end

        try
            if ~exist(dlgAns{1}, 'file')
                fclose(fopen(dlgAns{1}, 'w'));
            end
            hCfg_.setConfigFile(dlgAns{1}, 0);
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.rawRecordings = cellfun(@strip, strsplit(dlgAns{2}, ','), 'UniformOutput', 0);
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.sampleRate = str2double(dlgAns{3});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.nChans = str2double(dlgAns{4});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.bitScaling = str2double(dlgAns{5});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.headerOffset = str2double(dlgAns{6});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg_.dataType = dlgAns{7};
        catch ME
            errordlg(ME.message);
            continue;
        end

        break;
    end

    dlgAns = questdlg('Would you like to export advanced parameters as well?', 'Bootstrap', 'No');
    switch dlgAns
        case 'Yes'
            hCfg_.save(hCfg_.configFile, 1, 0);

        case 'No'
            hCfg_.save(hCfg_.configFile, 0, 0);

        otherwise
            return;
    end

    obj.hCfg = hCfg_;
    obj.hCfg.edit();
end

%% LOCAL FUNCTIONS
% function hCfg = bootstrapGUI() % WIP
%     %BOOTSTRAPGUI Show all (common) parameters
%     % load old2new param set and convert to new2old
%     [old2new, new2old] = jrclust.utils.getOldParamMapping();
% 
%     % build the bootstrap GUI
%     hBootstrap = uicontainer();
%     hRecData = uipanel('Parent', hBootstrap, 'Title', 'Recording file', 'Position', [0, 0.75, 0.25, 0.25]);
%     hProbe = uipanel('Parent', hBootstrap, 'Title', 'Probe parameters', 'Position', [0.25, 0.75, 0.25, 0.25]);
% end