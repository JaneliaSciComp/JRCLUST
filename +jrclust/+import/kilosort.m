function [hCfg, res] = kilosort(loadPath)
    %KILOSORT Import a Kilosort session from NPY files
    [hCfg, res] = deal([]);

    if exist('readNPY', 'file') ~= 2
        warning('Please make sure you have npy-matlab installed (https://github.com/kwikteam/npy-matlab)');
        return;
    end

    loadPath_ = jrclust.utils.absPath(loadPath);
    if isempty(loadPath_)
        error('Could not find path ''%s''', loadPath);
    elseif exist(loadPath, 'dir') ~= 7
        error('''%s'' is not a directory', loadPath);
    end

    loadPath = loadPath_;

    cfgData = struct();
    cfgData.outputDir = loadPath;
    
    % load params and set them in cfgData
    params = parseParams(fullfile(loadPath, 'params.py'));
    channelMap = readNPY(fullfile(loadPath, 'channel_map.npy')) + 1;
    channelPositions = readNPY(fullfile(loadPath, 'channel_positions.npy'));

    cfgData.sampleRate = params.sample_rate;
    cfgData.nChans = params.n_channels_dat;
    cfgData.dataType = params.dtype;
    cfgData.headerOffset = params.offset;
    cfgData.siteMap = channelMap;
    cfgData.siteLoc = channelPositions;
    cfgData.shankMap = ones(size(channelMap), 'like', channelMap); % this can change with a prm file
    cfgData.rawRecordings = {params.dat_path};

    hCfg = jrclust.Config(cfgData);
    
    hCfg.updateLog('import-kilosort', sprintf('Loading NPY files from %s', loadPath), 1, 0);

    % load spike data
    amplitudes = readNPY(fullfile(loadPath, 'amplitudes.npy'));
    spikeTimes = readNPY(fullfile(loadPath, 'spike_times.npy')) + 1;
    spikeTemplates = readNPY(fullfile(loadPath, 'spike_templates.npy')) + 1;
    spikeClusters = readNPY(fullfile(loadPath, 'spike_clusters.npy')) + 1;
    simScore = readNPY(fullfile(loadPath, 'similar_templates.npy'));
    templates = readNPY(fullfile(loadPath, 'templates.npy')); % nTemplates x nSamples x nChannels
    
    cProj = readNPY('template_features.npy')';
    iNeigh = readNPY('template_feature_ind.npy')';
    cProjPC = permute(readNPY('pc_features.npy'), [2 3 1]); % nFeatures x nSites x nSpikes
    iNeighPC = readNPY('pc_feature_ind.npy')';

    hCfg.updateLog('import-kilosort', 'Finished loading files', 0, 1);

    [clusterIDs, ~, indices] = unique(spikeClusters);
    goodClusters = clusterIDs(clusterIDs > 0);
    junkClusters = setdiff(clusterIDs, goodClusters);
    clusterIDsNew = [junkClusters' 1:numel(goodClusters)]';
    spikeClusters = clusterIDsNew(indices);

    nTemplates = size(templates, 1);
    nClusters = numel(goodClusters);

    spikeSites = zeros(size(spikeClusters), 'like', spikeClusters);
    for iTemplate = 1:nTemplates
        template = squeeze(templates(iTemplate, :, :));
        [~, tSite] = min(min(template));

        spikeSites(spikeTemplates == iTemplate) = tSite;
    end

    %%% try to detect the recording file
    % first check for a .meta file
    binfile = params.dat_path;
    metafile = jrclust.utils.absPath(jrclust.utils.subsExt(binfile, '.meta'));
    if isempty(metafile)
        dlgAns = questdlg('Do you have a .meta file?', 'Import', 'No');

        switch dlgAns
            case 'Yes' % select .meta file
                [metafile, loadPath] = jrclust.utils.selectFile({'*.meta', 'SpikeGLX meta files (*.meta)'; '*.*', 'All Files (*.*)'}, 'Select a .meta file', loadPath, 0);
                if isempty(metafile)
                    return;
                end

                if isempty(binfile)
                    binfile = jrclust.utils.subsExt(metafile, '.bin');
                end

            case 'No' % select recording file
                if isempty(binfile)
                    [binfile, loadPath] = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select a raw recording', loadPath, 0);
                    if isempty(binfile)
                        return;
                    end
                end

            case {'Cancel', ''}
                return;
        end
    end

    % check for missing binary file
    binfile = jrclust.utils.absPath(binfile);
    if isempty(jrclust.utils.absPath(binfile))
        binfile = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select a raw recording', loadPath, 0);
        if isempty(binfile)
            return;
        end
    end

    % load metafile, set bitScaling
    if ~isempty(metafile)
        SMeta_ = jrclust.utils.loadMetadata(metafile);
        hCfg.bitScaling = SMeta_.bitScaling;
    else
        hCfg.bitScaling = 1;
    end

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

    % set some specific params
    hCfg.nPeaksFeatures = 1; % don't find secondary peaks
    hCfg.figList = setdiff(hCfg.figList, 'FigRD'); % don't show rho-delta plot
    hCfg.corrRange = [0.75 1];

    %%% detect and extract spikes/features
    hDetect = jrclust.detect.DetectController(hCfg, spikeTimes, spikeSites);
    dRes = hDetect.detect();
    dRes.spikeSites = spikeSites;
    sRes = struct('spikeClusters', spikeClusters, ...
                  'spikeTemplates', spikeTemplates, ...
                  'simScore', simScore, ...
                  'amplitudes', amplitudes, ...
                  'templateFeatures', cProj, ...
                  'templateFeatureInd', iNeigh, ...
                  'pcFeatures', cProjPC, ...
                  'pcFeatureInd', iNeighPC);

    hClust = jrclust.sort.TemplateClustering(sRes, dRes, hCfg);
    hClust.computeCentroids();
    hClust.updateWaveforms();
    hClust.computeQualityScores();

    res = jrclust.utils.mergeStructs(dRes, sRes);
    res.hClust = hClust;
end

