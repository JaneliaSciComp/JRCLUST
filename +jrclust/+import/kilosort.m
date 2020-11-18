function [hCfg, res] = kilosort(loadPath,confirm_flag)
%KILOSORT Import a Kilosort session from NPY files

if nargin<2
    confirm_flag=true;
end
[hCfg, res] = deal([]);

phyData = loadPhy(loadPath);
if ~isfield(phyData, 'loadPath')
    return;
end
loadPath = phyData.loadPath;

cfgData = struct();

% load params and set them in cfgData
params = phyData.params;
channelMap = phyData.channel_map + 1;
channelPositions = phyData.channel_positions;

% check for existence of .prm file in one of two places:
% 1. loadPath (where the .npy files are found)
% 2. datPath (where the recording is found)
names = dir(loadPath);
loadPathPrmCandidate = contains({names.name}, '.prm');

[dirname, filename, ~] = fileparts(params.dat_path);
datPathPrmCandidate = fullfile(dirname, [filename, '.prm']);

if sum(loadPathPrmCandidate) == 1
    cfgData.template_file = fullfile(loadPath, names(loadPathPrmCandidate).name);
    cfgData = jrclust.utils.mergeStructs(cfgData, jrclust.utils.mToStruct(cfgData.templateFile));
elseif exist(datPathPrmCandidate, 'file') == 2
    cfgData.template_file = datPathPrmCandidate;
    cfgData = jrclust.utils.mergeStructs(cfgData, jrclust.utils.mToStruct(cfgData.templateFile));
else
    warning('Could not find a .prm file to use as a template.');
end

% fill in any gaps from Phy params
if ~isfield(cfgData, 'outputDir')
    cfgData.outputDir = loadPath;
end
if ~isfield(cfgData, 'sampleRate')
    cfgData.sampleRate = params.sample_rate;
end
if ~isfield(cfgData, 'nChans')
    cfgData.nChans = params.n_channels_dat;
end
if ~isfield(cfgData, 'dataTypeRaw')
    cfgData.dataTypeRaw = params.dtype;
end
if ~isfield(cfgData, 'dataTypeExtracted')
    cfgData.dataTypeExtracted = params.dtype;
end
if ~isfield(cfgData, 'headerOffset')
    cfgData.headerOffset = params.offset;
end
if ~isfield(cfgData, 'siteMap')
    cfgData.siteMap = channelMap;
end
if ~isfield(cfgData, 'channelPositions')
    cfgData.siteLoc = channelPositions;
end
if ~isfield(cfgData, 'shankMap')
    cfgData.shankMap = ones(size(channelMap), 'like', channelMap);
end
if ~isfield(cfgData, 'rawRecordings')
    cfgData.rawRecordings = {params.dat_path};
end

hCfg = jrclust.Config(cfgData);

siteMapIndex = find(ismember(hCfg.siteMap,channelMap)); % KS may have returned a partial channel map, i.e. ~all(ismember(hCfg.siteMap,channelMap)) 

% load spike data
amplitudes = phyData.amplitudes;
spikeTimes = phyData.spike_times + 1;
spikeTemplates = phyData.spike_templates + 1;
spikeClusters = phyData.spike_clusters + 1;
simScore = phyData.similar_templates;
templates = phyData.templates; % nTemplates x nSamples x nChannels

cProj = phyData.template_features';
iNeigh = phyData.template_feature_ind';
cProjPC = permute(phyData.pc_features, [2 3 1]); % nFeatures x nSites x nSpikes
iNeighPC = phyData.pc_feature_ind';

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
    [~, tSite] = max(max(abs(template)));
    spikeSites(spikeTemplates == iTemplate) = siteMapIndex(tSite); %setting spikeSites as tSites directly would renumber channel map in case KS has returned a partial channel map 
end

%%% try to detect the recording file
% first check for a .meta file
binfile = hCfg.rawRecordings{1};
if isempty(hCfg.rawRecordings{1})
    binfile = params.dat_path; % take ap.bin file location from .prm file by default
end
metafile = jrclust.utils.absPath(jrclust.utils.subsExt(binfile, '.meta'));
if isempty(metafile)
    dlgAns = questdlg('Do you have a .meta file?', 'Import', 'Yes','No','Cancel','Yes');

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
    % confirm with the user if confirm_flag is true
    [~, sessionName, ~] = fileparts(hCfg.rawRecordings{1});
    configFile = fullfile(hCfg.outputDir, [sessionName, '.prm']);

    dlgFieldNames = {'Config filename', ...
                     'Raw recording file(s)', ...
                     'Sampling rate (Hz)', ...
                     'Number of channels in file', ...
                     sprintf('%sV/bit', char(956)), ...
                     'Header offset (bytes)', ...
                     'Data Type in Raw file (int16, uint16, single, double)', ...
                     'Data Type Extracted (int16, uint16, single, double)'};
    dlgFieldVals = {configFile, ...
                    strjoin(hCfg.rawRecordings, ','), ...
                    num2str(hCfg.sampleRate), ...
                    num2str(hCfg.nChans), ...
                    num2str(hCfg.bitScaling), ...
                    num2str(hCfg.headerOffset), ...
                    hCfg.dataTypeRaw,...
                    hCfg.dataTypeExtracted};
    if confirm_flag
        dlgAns = inputdlg(dlgFieldNames, 'Does this look correct?', 1, dlgFieldVals, struct('Resize', 'on', 'Interpreter', 'tex'));
    else
        dlgAns = dlgFieldVals;
    end
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
        hCfg.dataTypeRaw = dlgAns{7};
    catch ME
        errordlg(ME.message);
        continue;
    end

    try
        hCfg.dataTypeExtracted = dlgAns{8};
    catch ME
        errordlg(ME.message);
        continue;
    end
 
    break;
end
    
% remove out-of-bounds spike times
d = dir(hCfg.rawRecordings{1});
nSamples = d.bytes / jrclust.utils.typeBytes(hCfg.dataTypeRaw) / hCfg.nChans;
oob = spikeTimes > nSamples;
if any(oob)
    warning('Removing %d/%d spikes after the end of the recording', sum(oob), numel(oob));
    spikeTimes = spikeTimes(~oob);
    spikeTemplates = spikeTemplates(~oob);
    spikeSites = spikeSites(~oob);
    spikeClusters = spikeClusters(~oob);
    cProj = cProj(:, ~oob);
    cProjPC = cProjPC(:, :, ~oob);
end

% set some specific params
hCfg.nPeaksFeatures = 1; % don't find secondary peaks
% remove FigRD
if ismember(hCfg.figList,'FigRD')
    keepFigIdx = ~ismember(hCfg.figList,'FigRD');
    hCfg.figList = hCfg.figList(keepFigIdx);
    if ~isempty(hCfg.figPos)
        hCfg.figPos = hCfg.figPos(keepFigIdx);
    end
end
hCfg.corrRange = [0.75 1];

% save out param file
hCfg.save();

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

hClust = jrclust.sort.TemplateClustering(hCfg, sRes, dRes);

res = jrclust.utils.mergeStructs(dRes, sRes);
res.hClust = hClust;
end  %fun
