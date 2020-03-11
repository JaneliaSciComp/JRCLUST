function bootstrap(obj, varargin)
    %BOOTSTRAP Bootstrap a JRCLUST session
    %   metafile: optional string; path (or glob) to meta file(s)
    if nargin > 1
        metafile_ = jrclust.utils.absPath(varargin{1});
        if isempty(metafile_) % TODO: warn?
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
            [~, ~, exts] = cellfun(@(f) fileparts(f), metafile, 'UniformOutput', 0);
            uniqueExts = unique(exts);
            if numel(uniqueExts) > 1
                error('Specify only a single file type');
            end

            ext = uniqueExts{:};
            switch lower(ext)
                case '.rhd'
                    obj.bootstrapIntan(metafile);
                    return;

                case '.meta'
                    binfile = cellfun(@(f) jrclust.utils.subsExt(f, '.bin'), metafile, 'UniformOutput', 0);
            end
        end
        % set whether to ask user input 
        if any(cellfun(@(x) strcmp(x,'-noconfirm'), varargin))
            ask=false;
        else
            ask=true;
        end
        % check whether user requires advanced parameters
        if any(cellfun(@(x) strcmp(x,'-advanced'), varargin))
            advanced=true;
        else
            advanced=false;
        end
    else
        metafile = '';
        workingdir = pwd();
    end

    % first check for a .meta file
    if isempty(metafile)
        [metafile, binfile, workingdir] = getMetafile(workingdir);
    end

    % check for missing binary files
    if any(cellfun(@(f) isempty(jrclust.utils.absPath(f)), binfile))
        binfile = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; ...
                                            '*.*', 'All Files (*.*)'}, 'Select one or more raw recordings', workingdir, 1);
        if cellfun(@isempty, binfile)
            return;
        end
    end

    if ~isempty(metafile) % load metafile
        cfgData = metaToConfig(metafile, binfile, workingdir);
    else % ask for a probe file instead
        cfgData = struct('outputDir', workingdir);
        cfgData.rawRecordings = binfile;
        cfgData.probe_file = getProbeFile(cfgData.outputDir,ask);

        if isempty(cfgData.probe_file) % closed dialog, cancel bootstrap
            return;
        end
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
                        num2str(hCfg_.sampleRate, 15), ...
                        num2str(hCfg_.nChans), ...
                        num2str(hCfg_.bitScaling), ...
                        num2str(hCfg_.headerOffset), ...
                        hCfg_.dataType};
        if ask
            dlgAns = inputdlg(dlgFieldNames, 'Does this look correct?', 1, dlgFieldVals, struct('Resize', 'on', 'Interpreter', 'tex'));
        else
            dlgAns= dlgFieldVals';
        end
        if isempty(dlgAns)
            return;
        end

        try
            if ~exist(dlgAns{1}, 'file')
                fclose(fopen(dlgAns{1}, 'w'));
            end
            hCfg_.setConfigFile(dlgAns{1}, 0);
            hCfg_.outputDir = fileparts(dlgAns{1}); % set outputdir to wherever configFile lives
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

    if ask
        dlgAns = questdlg('Would you like to export advanced parameters?', 'Bootstrap', 'No');
    elseif advanced
        dlgAns = 'Yes';
    else
        dlgAns = 'No';
    end
    
    switch dlgAns
        case 'Yes'
            hCfg_.save('', 1);

        case 'No'
            hCfg_.save('', 0);

        otherwise
            return;
    end

    obj.hCfg = hCfg_;
    obj.hCfg.edit();
end

%% LOCAL FUNCTIONS
function [metafile, binfile, workingdir] = getMetafile(workingdir)
    dlgAns = questdlg('Do you have a .meta file?', 'Bootstrap', 'No');
    switch dlgAns
        case 'Yes' % select .meta file
            [metafile, workingdir] = jrclust.utils.selectFile({'*.meta', 'SpikeGLX meta files (*.meta)'; '*.*', 'All Files (*.*)'}, 'Select one or more .meta files', workingdir, 1);
            if all(cellfun(@isempty, metafile))
                return;
            end

            binfile = cellfun(@(f) jrclust.utils.subsExt(f, '.bin'), metafile, 'UniformOutput', 0);

        case 'No' % select recording file
            metafile = '';
            [binfile, workingdir] = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; ...
                                                              '*.rhd', 'Intan recordings (*.rhd)'; ...
                                                              '*.*', 'All Files (*.*)'}, 'Select one or more raw recordings', workingdir, 1);
            if all(cellfun(@isempty, binfile))
                return;
            end

            [~, ~, exts] = cellfun(@(f) fileparts(f), binfile, 'UniformOutput', 0);
            uniqueExts = unique(exts);
            if numel(uniqueExts) > 1
                error('Specify only a single file type');
            end

            ext = uniqueExts{:};
            if strcmpi(ext, '.rhd')
                obj.bootstrapIntan(binfile);
                return;
            end

        case {'Cancel', ''}
            return;
    end
end

function cfgData = metaToConfig(metafile, binfile, workingdir)
    cfgData = struct('outputDir', workingdir);
    cfgData.rawRecordings = binfile;

    SMeta = jrclust.utils.loadMetadata(metafile{1});
    cfgData.sampleRate = SMeta.sampleRate;
    cfgData.nChans = SMeta.nChans;
    cfgData.bitScaling = SMeta.bitScaling;
    cfgData.headerOffset = 0; % standard for SpikeGLX
    cfgData.dataType = 'int16'; % standard for SpikeGLX

    probeFile = getProbeFile(cfgData.outputDir, 0);

    if isempty(probeFile) && SMeta.isImec % think we've got a Neuropixels probe
        if ~isempty(SMeta.probeOpt) % 3A with option
            probeFile = sprintf('imec3_opt%d.prb', SMeta.probeOpt);
        else % 3A or 3B; ask
            dlgAns = questdlg('It looks like you have a Neuropixels probe. Is this correct?', 'Bootstrap', ...
                              'Yes', 'No', 'No'); % don't permit 'Cancel'
            if isempty(dlgAns) % closed dialog; cancel
                cfgData = [];
                return;
            end

            if strcmp(dlgAns, 'Yes')
                dlgAns = listdlg('PromptString', 'Specify your configuration', ...
                                 'SelectionMode', 'single', ...
                                 'ListString', {'Phase 3A', ...
                                                'Phase 3B (Staggered)', ...
                                                'Phase 3B (Aligned)', ...
                                                'Custom configuration'});
                switch dlgAns
                    case 1
                        probeFile = 'imec3a.prb';
                    case 2
                        probeFile = 'imec3b_staggered.prb';
                    case 3
                        probeFile = 'imec3b_aligned.prb';
                    case 4
                        probeFile = getProbeFile(cfgData.outputDir);
                end
            end
        end
    elseif isempty(probeFile)
        probeFile = getProbeFile(cfgData.outputDir);
    end

    cfgData.probe_file = probeFile;
end

function probeFile = probeFileInDir(workingdir)
    probeFile = '';

    d = dir(fullfile(workingdir, '*.prb'));
    if isempty(d)
        return;
    end

    probeFile = fullfile(workingdir, d(1).name);
end

function probeFile = getProbeFile(workingdir, ask)
    if nargin < 2
        ask = 1;
    end

    probeFile = probeFileInDir(workingdir);

    if ~isempty(probeFile) && ask % found a probe file in working directory; confirm
        dlgAns = questdlg(sprintf('Found probe file ''%s''. Use it?', probeFile));
        if strcmp(dlgAns, 'Yes')
            return;
        else
            probeFile = '';
        end
    end

    if ask
        dlgAns = questdlg('Would you like to specify a probe file?', 'Bootstrap', 'No');
        if strcmp(dlgAns, 'Yes')
            probedir = workingdir;
            if isempty(dir(fullfile(workingdir, '*.prb')))
                probedir = fullfile(jrclust.utils.basedir(), 'probes');
            end

            [probefile, probedir] = jrclust.utils.selectFile({'*.prb', 'Probe files (*.prb)'; '*.*', 'All Files (*.*)'}, ...
                                                             'Select a probe file', probedir, 0);
            probeFile = fullfile(probedir, probefile);
        end
    end
end
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