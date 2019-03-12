function bootstrapIntan(obj, rhdFiles)
    %BOOTSTRAPINTAN Summary of this function goes here
    cfgData = struct();

    for iFile = 1:numel(rhdFiles)
        rawPath = jrclust.utils.absPath(rhdFiles{iFile});
        if isempty(rawPath)
            continue;
        end

        fileData(iFile) = jrclust.detect.IntanRecording.readHeader(rawPath);
    end

    workingdir = fileparts(rhdFiles{1});

    % check some key fields
    nChans = {fileData.nChans};
    if numel(unique([nChans{:}])) == 1
        cfgData.nChans = nChans{1};
    else
        error('Recordings must have the same channel count');
    end

    sampleRate = {fileData.sampleRate};
    if numel(unique([sampleRate{:}])) == 1
        cfgData.sampleRate = sampleRate{1};
    else
        error('Recordings sampled at different rates');
    end

    headerOffset = {fileData.headerOffset};
    if numel(unique([headerOffset{:}])) == 1
        cfgData.headerOffset = headerOffset{1};
    else
        cfgData.headerOffset = 0;
    end

    cfgData.dataType = 'single';
    cfgData.bitScaling = 1;
    cfgData.rawRecordings = rhdFiles;
    cfgData.recordingFormat = 'Intan';

    dlgAns = questdlg('Would you like to specify a probe file?', 'Bootstrap', 'Yes');
    switch dlgAns
        case 'Yes' % select .prb file
            probedir = workingdir;
            if isempty(dir(fullfile(workingdir, '*.prb')))
                probedir = fullfile(jrclust.utils.basedir(), 'probes');
            end
            [probefile, probedir] = jrclust.utils.selectFile({'*.prb', 'Probe files (*.prb)'; '*.*', 'All Files (*.*)'}, 'Select a probe file', probedir, 0);
            cfgData.probe_file = fullfile(probedir, probefile);
            
        case 'No'
            % set some trivial defaults
            cfgData.siteMap = 1:cfgData.nChans;
            cfgData.siteLoc = [zeros(cfgData.nChans, 1), 50*(0:cfgData.nChans-1)'];
            cfgData.shankMap = ones(cfgData.nChans, 1);

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

    dlgAns = questdlg('Would you like to export advanced parameters as well?', 'Bootstrap', 'No');
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

