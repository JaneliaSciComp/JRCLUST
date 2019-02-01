function hCfg = doBootstrap(varargin)
    %DOBOOTSTRAP Bootstrap a JRCLUST session
    %   metafile: optional string; path (or glob) to meta file(s)
    if nargin > 0
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

    hCfg = [];

    % first check for a .meta file
    if isempty(metafile)
        dlgAns = questdlg('Do you have a .meta file?', 'Bootstrap', 'No');

        switch dlgAns
            case 'Yes' % select .meta file
                [metafile, workingdir] = uigetfile({'*.meta', 'SpikeGLX meta files (*.meta)'; '*.*', 'All Files (*.*)'}, 'Select one or more .meta files', workingdir, 'MultiSelect', 'on');
                if ~(ischar(metafile) || iscell(metafile)) % cancel
                    return;
                elseif ischar(metafile)
                    metafile = {fullfile(workingdir, metafile)};
                else
                    metafile = cellfun(@(f) fullfile(workingdir, f), metafile, 'UniformOutput', 0);
                end

                binfile = cellfun(@(f) jrclust.utils.subsExt(f, '.bin'), metafile, 'UniformOutput', 0);

            case 'No' % select recording file
                [binfile, workingdir] = uigetfile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select one or more raw recordings', workingdir, 'MultiSelect', 'on');
                if isempty(binfile) % cancel
                    return;
                end
                if ~(ischar(binfile) || iscell(binfile)) % cancel
                    return;
                elseif ischar(binfile)
                    binfile = {fullfile(workingdir, binfile)};
                else
                    binfile = cellfun(@(f) fullfile(workingdir, f), binfile, 'UniformOutput', 0);
                end

            case {'Cancel', ''}
                hCfg = [];
                return;
        end
    end

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
            [probefile, probedir] = uigetfile({'*.prb', 'Probe files (*.prb)'; '*.*', 'All Files (*.*)'}, 'Select a probe file', probedir);
            cfgData.probe_file = fullfile(probedir, probefile);

        case {'Cancel', ''}
            hCfg = [];
            return;
    end

    % construct the Config object from specified data
    hCfg = jrclust.Config(cfgData);

    while 1
        % confirm with the user
        [~, sessionName, ~] = fileparts(hCfg.rawRecordings{1});
        configFile = fullfile(hCfg.outputDir, [sessionName, '.prm']);

        dlgFieldNames = {'Config filename', ...
                         'Raw recording file(s)', ...
                         'Sampling rate (Hz)', ...
                         'Number of channels in file', ...
                         sprintf('%sV/bit', char(956)), ...
                         'Header offset (bytes)', ...
                         'Data Type (int16, uint16, single, double)'};
        dlgFieldVals = {configFile, ...
                        strjoin(hCfg.rawRecordings, ','), ...
                        num2str(hCfg.sampleRate), ...
                        num2str(hCfg.nChans), ...
                        num2str(hCfg.bitScaling), ...
                        num2str(hCfg.headerOffset), ...
                        hCfg.dataType};
        dlgAns = inputdlg(dlgFieldNames, 'Does this look correct?', 1, dlgFieldVals, struct('Resize', 'on', 'Interpreter', 'tex'));
        if isempty(dlgAns)
            return;
        end

        try
            if ~exist(dlgAns{1}, 'file')
                fclose(fopen(dlgAns{1}, 'w'));
            end
            hCfg.setConfigFile(dlgAns{1}, 0);
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.rawRecordings = cellfun(@strip, strsplit(dlgAns{2}, ','), 'UniformOutput', 0);
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.sampleRate = str2double(dlgAns{3});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.nChans = str2double(dlgAns{4});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.bitScaling = str2double(dlgAns{5});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.headerOffset = str2double(dlgAns{6});
        catch ME
            errordlg(ME.message);
            continue;
        end

        try
            hCfg.dataType = dlgAns{7};
        catch ME
            errordlg(ME.message);
            continue;
        end

        break;
    end

    dlgAns = questdlg('Would you like to export advanced parameters as well?', 'Bootstrap', 'No');
    switch dlgAns
        case 'Yes'
            hCfg.save(hCfg.configFile, 1, 0);

        case 'No'
            hCfg.save(hCfg.configFile, 0, 0);

        otherwise
            return;
    end

    hCfg.edit();
end

%% LOCAL FUNCTIONS
function hCfg = bootstrapGUI() % WIP
    %BOOTSTRAPGUI Show all (common) parameters
    % load old2new param set and convert to new2old
    [old2new, new2old] = jrclust.utils.getOldParamMapping();

    % build the bootstrap GUI
    hBootstrap = uicontainer();
    hRecData = uipanel('Parent', hBootstrap, 'Title', 'Recording file', 'Position', [0, 0.75, 0.25, 0.25]);
    hProbe = uipanel('Parent', hBootstrap, 'Title', 'Probe parameters', 'Position', [0.25, 0.75, 0.25, 0.25]);


end
