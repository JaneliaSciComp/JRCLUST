function [hCfg, res] = kilosort(rezFile)
    %KILOSORT Import a Kilosort session from rez.mat
    [hCfg, res] = deal([]);

    if exist('readNPY', 'file') ~= 2
        warning('Please make sure you have npy-matlab installed (https://github.com/kwikteam/npy-matlab)');
        return;
    end

    rezFile_ = jrclust.utils.absPath(rezFile);
    if isempty(rezFile_)
        error('Could not find file ''%s''', rezFile);
    end

    rezFile = rezFile_;
    workingdir = fileparts(rezFile);

    try
        load(rezFile, 'rez');
    catch ME
        error('Failed to load ''%s'': %s', rezFile, ME.message);
    end

    rezFields = fieldnames(rez); %#ok<NODEF>
    for i = 1:numel(rezFields)
        fn = rezFields{i};
        if isa(rez.(fn), 'gpuArray')
            rez.(fn) = jrclust.utils.tryGather(rez.(fn));
        end
    end

    spikeTimes = rez.st3(:, 1);
    spikeTemplates = rez.st3(:, 2);

    amplitudes = rez.st3(:, 3);

    if ~isfield(rez, 'ccb') && size(rez.st3, 2) > 4 % Kilosort1
        spikeClusters = 1 + rez.st3(:, 5);
    else
        spikeClusters = spikeTemplates;
    end

    [clusterIDs, ~, indices] = unique(spikeClusters);
    goodClusters = clusterIDs(clusterIDs > 0);
    junkClusters = setdiff(clusterIDs, goodClusters);
    clusterIDsNew = [junkClusters' 1:numel(goodClusters)]';
    spikeClusters = clusterIDsNew(indices);
    nTemplates = size(rez.simScore, 1);

    % nClusters = numel(goodClusters);

    % clusterTemplates = arrayfun(@(iCluster) unique(spikeTemplates(spikeClusters == iCluster)), 1:nClusters, 'UniformOutput', 0);

    % compute templates
    nt0 = size(rez.W, 1);
    U = rez.U;
    W = rez.W;

    Nfilt = size(W, 2);
    Nchan = rez.ops.Nchan;

    templates = zeros(Nchan, nt0, Nfilt, 'single');
    for iNN = 1:size(templates,3)
       templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
    end
    templates = permute(templates, [3 2 1]); % nTemplates x nSamples x nChannels

    spikeSites = zeros(size(spikeClusters), 'like', spikeClusters);
    for iTemplate = 1:nTemplates
        template = squeeze(templates(iTemplate, :, :));
        [~, tSite] = min(min(template));

        spikeSites(spikeTemplates == iTemplate) = tSite;
    end

    % import features
    if isfield(rez, 'cProj') && ~isempty(rez.cProj)
        cProj = rez.cProj;
        iNeigh = rez.iNeigh;
    elseif exist(fullfile(workingdir, 'template_features.npy'), 'file') == 2
        cProj = readNPY('template_features.npy')';
        iNeigh = readNPY('template_feature_ind.npy')';
    else
        [cProj, iNeigh] = deal([]);
    end

    if isfield(rez, 'cProjPC') && ~isempty(rez.cProjPC)
        cProjPC = permute(rez.cProjPC, [2 3 1]); % nFeatures x nSites x nSpikes
        iNeighPC = rez.iNeighPC;
    elseif exist(fullfile(workingdir, 'pc_features.npy'), 'file') == 2
        cProjPC = permute(readNPY('pc_features.npy'), [2 3 1]); % nFeatures x nSites x nSpikes
        iNeighPC = readNPY('pc_feature_ind.npy')';
    else
        [cProjPC, iNeighPC] = deal([]);
    end

    %%% try to detect the recording file
    % first check for a .meta file
    binfile = jrclust.utils.absPath(rez.ops.fbinary);
    metafile = jrclust.utils.absPath(jrclust.utils.subsExt(binfile, '.meta'));
    if isempty(metafile)
        dlgAns = questdlg('Do you have a .meta file?', 'Import', 'No');

        switch dlgAns
            case 'Yes' % select .meta file
                [metafile, workingdir] = jrclust.utils.selectFile({'*.meta', 'SpikeGLX meta files (*.meta)'; '*.*', 'All Files (*.*)'}, 'Select a .meta file', workingdir, 0);
                if isempty(metafile)
                    return;
                end

                if isempty(binfile)
                    binfile = jrclust.utils.subsExt(metafile, '.bin');
                end

            case 'No' % select recording file
                if isempty(binfile)
                    [binfile, workingdir] = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select a raw recording', workingdir, 0);
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
        binfile = jrclust.utils.selectFile({'*.bin;*.dat', 'SpikeGLX recordings (*.bin, *.dat)'; '*.*', 'All Files (*.*)'}, 'Select a raw recording', workingdir, 0);
        if isempty(binfile)
            return;
        end
    end

    % load metafile
    if ~isempty(metafile)
        SMeta_ = jrclust.utils.loadMetadata(metafile);
        cfgData = struct('sampleRate', SMeta_.sampleRate, ...
                         'nChans', SMeta_.nChans, ...
                         'bitScaling', SMeta_.bitScaling, ...
                         'headerOffset', 0, ...
                         'dataType', SMeta_.dataType);

        cfgData.rawRecordings = binfile;
        cfgData.outputDir = workingdir;
    else
        cfgData.sampleRate = rez.ops.fs;
        cfgData.nChans = rez.ops.NchanTOT;
        cfgData.dataType = 'int16';
        cfgData.rawRecordings = binfile;
        cfgData.outputDir = workingdir;
    end

    dlgAns = questdlg('Would you like to specify a probe file?', 'Import', 'No');
    switch dlgAns
        case 'Yes' % select .prb file
            probedir = workingdir;
            if isempty(dir(fullfile(workingdir, '*.prb')))
                probedir = fullfile(jrclust.utils.basedir(), 'probes');
            end
            [probefile, probedir] = jrclust.utils.selectFile({'*.prb', 'Probe files (*.prb)'; '*.*', 'All Files (*.*)'}, 'Select a probe file', probedir, 0);
            cfgData.probe_file = fullfile(probedir, probefile);

        case 'No'
            if isfield(rez, 'connected')
                cfgData.siteMap = rez.ops.chanMap(rez.connected);
                if isfield(rez, 'xc')
                    cfgData.siteLoc = [rez.xc(:) rez.yc(:)];
                else
                    cfgData.siteLoc = [rez.xcoords(:) rez.ycoords(:)];
                    cfgData.siteLoc = cfgData.siteLoc(connected, :);
                end
            else
                cfgData.siteMap = rez.ops.chanMap;
                cfgData.siteLoc = [rez.xc(:) rez.yc(:)];
            end

            cfgData.shankMap = rez.ops.kcoords;

        case {'Cancel', ''}
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
                  'simScore', rez.simScore, ...
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

