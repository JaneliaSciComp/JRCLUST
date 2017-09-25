function varargout = kilosort(varargin)
% Marius Pachitariu Kilosort
addpath('./kilosort');
addpath('./npy-matlab');

% kilosort
vcMode = varargin{1};
switch vcMode
    case 'config'
        varargout{1} = config_(varargin{2});
    case 'chanMap'
        varargout{1} = chanMap_(varargin{2}); %converts .prm file to chanMap format\
    case 'rezToPhy'
        rezToPhy_(varargin{2}, varargin{3});
    case 'preprocessData'
        [varargout{1}, varargout{2}, varargout{3}] = preprocessData_(varargin{2});
    case 'fitTemplates'
        varargout{1} = fitTemplates_(varargin{2}, varargin{3}, varargin{4});
    case 'fullMPMU'
        varargout{1} = fullMPMU_(varargin{2}, varargin{3});
    case 'merge_posthoc2'
        varargout{1} = merge_posthoc2_(varargin{2});
    case 'merge_posthoc2_jrc2'
        varargout{1} = merge_posthoc2_jrc2_(varargin{2});
end %switch        
end %func


%--------------------------------------------------------------------------
function ops = config_(P)
% ops = config_(P)
% ops = config_()
Nfilt_per_site = 2;
if nargin>=1
    useGPU = P.fGpu;
    [fpath, ~, ~] = fileparts(P.vcFile);
    nSites = numel(P.viSite2Chan);
    ops.NchanTOT = P.nChans; % total number of chanels stored infile
elseif nargin==0
    useGPU = 1;
    fpath = ''; %pwd
    nSites = nan;
end

ops.GPU                 = useGPU; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)		
ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm		
ops.verbose             = 1; % whether to print command line progress		
ops.showfigures         = 1; % whether to plot figures during optimization		
		
ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'		
ops.fbinary             = P.vcFile; %fullfile(fpath, 'sim_binary.dat'); % will be created for 'openEphys'		
ops.fproc               = fullfile(fpath, 'temp_wh.dat'); % residual from RAM of preprocessed data		
ops.root                = fpath; % 'openEphys' only: where raw files are		
% define the channel map as a filename (string) or simply an array		
ops.chanMap             = fullfile(fpath, 'chanMap.mat'); % make this file using createChannelMapFile.m		
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if unavailable chanMap file		

% ops.Nfilt               = ceil(nSites*Nfilt_per_site/32)*32;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.Nfilt               = ceil(nSites*Nfilt_per_site/32)*32;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)     		
ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)		
ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)		
		
% options for channel whitening		
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)		
ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)		
ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)		
		
ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). 		

% other options for controlling the model and optimization		
ops.Nrank               = 3;    % matrix rank of spike template model (3)		
ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)		
ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)		
ops.fshigh              = 200;   % frequency for high pass filtering		
% ops.fslow             = 2000;   % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection		
ops.scaleproc           = 200;   % int16 scaling of whitened data		
ops.NT                  = 128*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory) 		
% for GPU should be multiple of 32 + ntbuff		
		
% the following options can improve/deteriorate results. 		
% when multiple values are provided for an option, the first two are beginning and ending anneal values, 		
% the third is the value used in the final pass. 		
ops.Th               = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])		
ops.lam              = [5 5 5];   % large means amplitudes are forced around the mean ([10 30 30])		
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)		
ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])		
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)		
ops.mergeT           = .1;           % upper threshold for merging (.1)		
ops.splitT           = .1;           % lower threshold for splitting (.1)		
		
% options for initializing spikes from data		
ops.initialize      = 'fromData';    %'fromData' or 'no'		
ops.spkTh           = -6;      % spike threshold in standard deviations (4)		
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])		
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])		
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])		
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)		
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)		
		
% load predefined principal components (visualization only (Phy): used for features)		
dd                  = load('./kilosort/PCspikes2.mat'); % you might want to recompute this from your own data		
ops.wPCA            = dd.Wi(:,1:7);   % PCs 		
		
% options for posthoc merges (under construction)		
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)		
ops.epu     = Inf;		
		
ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
end %func


%--------------------------------------------------------------------------
function S_chanMap = chanMap_(P)
[fpath, ~, ~] = fileparts(P.vcFile);
chanMap = P.viSite2Chan;
chanMap = chanMap(:)';
% chanMap = [33 34 8 10 12 14 16 18 20 22 24 26 28 30 32 ...
%     7 9 11 13 15 17 19 21 23 25 27 29 31 1 2 3 4 5 6];

% the first thing Kilosort does is reorder the data with data = data(chanMap, :).
% Now we declare which channels are "connected" in this normal ordering, 
% meaning not dead or used for non-ephys data

% connected = true(34, 1); connected(1:2) = 0;
connected = true(size(chanMap));
connected(P.viSiteZero) = 0;

% now we define the horizontal (x) and vertical (y) coordinates of these
% 34 channels. For dead or nonephys channels the values won't matter. Again
% I will take this information from the specifications of the probe. These
% are in um here, but the absolute scaling doesn't really matter in the
% algorithm. 

% xcoords = 20 * [NaN NaN  1 0 0 1 0 1 0 1 0 1 0 1 0 1 0 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0];
% ycoords = 20 * [NaN NaN  7 8 9 9 10 10 11 11 12 12 13 13 14 14 15 15 16 ...
%     17 17 18 18 19 19 20 20 21 21 22 22 23 23 24]; 
xcoords = P.mrSiteXY(:,1)';
ycoords = P.mrSiteXY(:,2)';

% Often, multi-shank probes or tetrodes will be organized into groups of
% channels that cannot possibly share spikes with the rest of the probe. This helps
% the algorithm discard noisy templates shared across groups. In
% this case, we set kcoords to indicate which group the channel belongs to.
% In our case all channels are on the same shank in a single group so we
% assign them all to group 1. 

% kcoords = [NaN NaN 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
if isfield(P, 'shank')
    kcoords = P.shank;
else
    kcoords = ones(size(chanMap));
end
kcoords(P.viSiteZero) = nan;

% at this point in Kilosort we do data = data(connected, :), ycoords =
% ycoords(connected), xcoords = xcoords(connected) and kcoords =
% kcoords(connected) and no more channel map information is needed (in particular
% no "adjacency graphs" like in KlustaKwik). 
% Now we can save our channel map for the eMouse. 

% would be good to also save the sampling frequency here
% fs = 25000; 
fs = P.sRateHz;

S_chanMap = struct('chanMap', chanMap, 'connected', connected, 'xcoords', xcoords, 'ycoords', ycoords, 'kcoords', kcoords, 'fs', fs);

struct_save_(S_chanMap, fullfile(fpath, 'chanMap.mat'));

end %func


%--------------------------------------------------------------------------
function [spikeTimes, clusterIDs, amplitudes, templates, templateFeatures, ...
    templateFeatureInds, pcFeatures, pcFeatureInds] = rezToPhy_(rez, savePath)
% pull out results from kilosort's rez to either return to workspace or to
% save in the appropriate format for the phy GUI to run on. If you provide
% a savePath it should be a folder, and you will need to have npy-matlab
% available (https://github.com/kwikteam/npy-matlab)
%
% spikeTimes will be in samples, not seconds

fs = dir(fullfile(savePath, '*.npy'));
for i = 1:length(fs)
   delete(fullfile(savePath, fs(i).name)); 
end
if exist(fullfile(savePath, '.phy'), 'dir')
    rmdir(fullfile(savePath, '.phy'), 's');
end

spikeTimes = uint64(rez.st3(:,1));
% [spikeTimes, ii] = sort(spikeTimes);
spikeTemplates = uint32(rez.st3(:,2));
if size(rez.st3,2)>4
    spikeClusters = uint32(1+rez.st3(:,5));
end
amplitudes = rez.st3(:,3);

Nchan = rez.ops.Nchan;

% try
%     load(rez.ops.chanMap);
% catch
%    chanMap0ind  = [0:Nchan-1]';
%    connected    = ones(Nchan, 1);
%    xcoords      = ones(Nchan, 1);
%    ycoords      = (1:Nchan)';
% end
% chanMap0 = chanMap(connected>1e-6);

connected   = rez.connected(:);
xcoords     = rez.xcoords(:);
ycoords     = rez.ycoords(:);
chanMap     = rez.ops.chanMap(:);
chanMap0ind = chanMap - 1;

nt0 = size(rez.W,1);
U = rez.U;
W = rez.W;

% for i = 1:length(chanMap0)
%     chanMap0(i) = chanMap0(i) - sum(chanMap0(i) > chanMap(connected<1e-6));
% end
% [~, invchanMap0] = sort(chanMap0);

templates = zeros(Nchan, nt0, rez.ops.Nfilt, 'single');
for iNN = 1:rez.ops.Nfilt
   templates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))'; 
end
templates = permute(templates, [3 2 1]); % now it's nTemplates x nSamples x nChannels
templatesInds = repmat([0:size(templates,3)-1], size(templates,1), 1); % we include all channels so this is trivial

templateFeatures = rez.cProj;
templateFeatureInds = uint32(rez.iNeigh);
pcFeatures = rez.cProjPC;
pcFeatureInds = uint32(rez.iNeighPC);

if ~isempty(savePath)    
    writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
    writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_templates.npy')); % -1 for zero indexing
    if size(rez.st3,2)>4
        writeNPY(uint32(spikeClusters-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    else
        writeNPY(uint32(spikeTemplates-1), fullfile(savePath, 'spike_clusters.npy')); % -1 for zero indexing
    end
    writeNPY(amplitudes, fullfile(savePath, 'amplitudes.npy'));
    writeNPY(templates, fullfile(savePath, 'templates.npy'));
    writeNPY(templatesInds, fullfile(savePath, 'templates_ind.npy'));
    
%     Fs = rez.ops.fs;
    conn        = logical(connected);
    chanMap0ind = int32(chanMap0ind);
    
    writeNPY(chanMap0ind(conn), fullfile(savePath, 'channel_map.npy'));
    %writeNPY(connected, fullfile(savePath, 'connected.npy'));
%     writeNPY(Fs, fullfile(savePath, 'Fs.npy'));
    writeNPY([xcoords(conn) ycoords(conn)], fullfile(savePath, 'channel_positions.npy'));
    
    writeNPY(templateFeatures, fullfile(savePath, 'template_features.npy'));
    writeNPY(templateFeatureInds'-1, fullfile(savePath, 'template_feature_ind.npy'));% -1 for zero indexing
    writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
    writeNPY(pcFeatureInds'-1, fullfile(savePath, 'pc_feature_ind.npy'));% -1 for zero indexing
    
    whiteningMatrix = rez.Wrot/200;
    whiteningMatrixInv = whiteningMatrix^-1;
    writeNPY(whiteningMatrix, fullfile(savePath, 'whitening_mat.npy'));
    writeNPY(whiteningMatrixInv, fullfile(savePath, 'whitening_mat_inv.npy'));
    
    if isfield(rez, 'simScore')
        similarTemplates = rez.simScore;
        writeNPY(similarTemplates, fullfile(savePath, 'similar_templates.npy'));
    end
    
     %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        
        [~, fname, ext] = fileparts(rez.ops.fbinary);
        
        fprintf(fid,['dat_path = ''',fname ext '''\n']);
        fprintf(fid,'n_channels_dat = %i\n',rez.ops.NchanTOT);
        fprintf(fid,'dtype = ''int16''\n');
        fprintf(fid,'offset = 0\n');
        if mod(rez.ops.fs,1)
            fprintf(fid,'sample_rate = %i\n',rez.ops.fs);
        else
            fprintf(fid,'sample_rate = %i.\n',rez.ops.fs);
        end
        fprintf(fid,'hp_filtered = False');
        fclose(fid);
    end
end
end %func


%--------------------------------------------------------------------------
function [rez, DATA, uproj] = preprocessData_(ops)
tic;
uproj = [];
ops.nt0 	= getOr_(ops, {'nt0'}, 61);


if strcmp(ops.datatype , 'openEphys')
   ops = convertOpenEphysToRawBInary_(ops);  % convert data, only for OpenEphys
end

if ~isempty(ops.chanMap)
    if ischar(ops.chanMap)
        load(ops.chanMap);
        try
            chanMapConn = chanMap(connected>1e-6);
            xc = xcoords(connected>1e-6);
            yc = ycoords(connected>1e-6);
        catch
            chanMapConn = 1+chanNums(connected>1e-6);
            xc = zeros(numel(chanMapConn), 1);
            yc = [1:1:numel(chanMapConn)]';
        end
        ops.Nchan    = getOr_(ops, 'Nchan', sum(connected>1e-6));
        ops.NchanTOT = getOr_(ops, 'NchanTOT', numel(connected));
        if exist('fs', 'var')
            ops.fs       = getOr_(ops, 'fs', fs);
        end
    else
        chanMap = ops.chanMap;
        chanMapConn = ops.chanMap;
        xc = zeros(numel(chanMapConn), 1);
        yc = [1:1:numel(chanMapConn)]';
        connected = true(numel(chanMap), 1);      
        
        ops.Nchan    = numel(connected);
        ops.NchanTOT = numel(connected);
    end
else
    chanMap  = 1:ops.Nchan;
    connected = true(numel(chanMap), 1);
    
    chanMapConn = 1:ops.Nchan;    
    xc = zeros(numel(chanMapConn), 1);
    yc = [1:1:numel(chanMapConn)]';
end
if exist('kcoords', 'var')
    kcoords = kcoords(connected);
else
    kcoords = ones(ops.Nchan, 1);
end
NchanTOT = ops.NchanTOT;
NT       = ops.NT ;

rez.ops         = ops;
rez.xc = xc;
rez.yc = yc;
if exist('xcoords')
   rez.xcoords = xcoords;
   rez.ycoords = ycoords;
else
   rez.xcoords = xc;
   rez.ycoords = yc;
end
rez.connected   = connected;
rez.ops.chanMap = chanMap;
rez.ops.kcoords = kcoords; 

d = dir(ops.fbinary);
ops.sampsToRead = floor(d.bytes/NchanTOT/2);

if ispc
    dmem         = memory;
    memfree      = dmem.MemAvailableAllArrays/8;
    memallocated = min(ops.ForceMaxRAMforDat, dmem.MemAvailableAllArrays) - memfree;
    memallocated = max(0, memallocated);
else
    memallocated = ops.ForceMaxRAMforDat;
end
nint16s      = memallocated/2;

NTbuff      = NT + 4*ops.ntbuff;
Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));
Nbatch_buff = floor(4/5 * nint16s/rez.ops.Nchan /(NT-ops.ntbuff)); % factor of 4/5 for storing PCs of spikes
Nbatch_buff = min(Nbatch_buff, Nbatch);


%-----
% load data into patches, filter, compute covariance
if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
    [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
else
    [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
end

fprintf('Time %3.0fs. Loading raw data... \n', toc);
fid = fopen(ops.fbinary, 'r');
ibatch = 0;
Nchan = rez.ops.Nchan;
if ops.GPU
    CC = gpuArray.zeros( Nchan,  Nchan, 'single');
else
    CC = zeros( Nchan,  Nchan, 'single');
end
if strcmp(ops.whitening, 'noSpikes')
    if ops.GPU
        nPairs = gpuArray.zeros( Nchan,  Nchan, 'single');
    else
        nPairs = zeros( Nchan,  Nchan, 'single');
    end
end
if ~exist('DATA', 'var')
    DATA = zeros(NT, rez.ops.Nchan, Nbatch_buff, 'int16');
end

isproc = zeros(Nbatch, 1);
while 1
    ibatch = ibatch + ops.nSkipCov;
    
    offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
    if ibatch==1
        ioffset = 0;
    else
        ioffset = ops.ntbuff;
    end
    fseek(fid, offset, 'bof');
    buff = fread(fid, [NchanTOT NTbuff], '*int16');
    
    %         keyboard;
    
    if isempty(buff)
        break;
    end
    nsampcurr = size(buff,2);
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
    end
    if ops.GPU
        dataRAW = gpuArray(buff);
    else
        dataRAW = buff;
    end
    dataRAW = dataRAW';
    dataRAW = single(dataRAW);
    dataRAW = dataRAW(:, chanMapConn);
    
    datr = filter(b1, a1, dataRAW);
    datr = flipud(datr);
    datr = filter(b1, a1, datr);
    datr = flipud(datr);
    
    switch ops.whitening
        case 'noSpikes'
            smin      = my_min_(datr, ops.loc_range, [1 2]);
            sd = std(datr, [], 1);
            peaks     = single(datr<smin+1e-3 & bsxfun(@lt, datr, ops.spkTh * sd));
            blankout  = 1+my_min_(-peaks, ops.long_range, [1 2]);
            smin      = datr .* blankout;
            CC        = CC + (smin' * smin)/NT;
            nPairs    = nPairs + (blankout'*blankout)/NT;
        otherwise
            CC        = CC + (datr' * datr)/NT;
    end
    
    if ibatch<=Nbatch_buff
        DATA(:,:,ibatch) = gather_try_(int16( datr(ioffset + (1:NT),:)));
        isproc(ibatch) = 1;
    end
end
CC = CC / ceil((Nbatch-1)/ops.nSkipCov);
switch ops.whitening
    case 'noSpikes'
        nPairs = nPairs/ibatch;
end
fclose(fid);
fprintf('Time %3.0fs. Channel-whitening filters computed. \n', toc);
switch ops.whitening
    case 'diag'
        CC = diag(diag(CC));
    case 'noSpikes'
        CC = CC ./nPairs;
end

if ops.whiteningRange<Inf
    ops.whiteningRange = min(ops.whiteningRange, Nchan);
    Wrot = whiteningLocal_(gather_try_(CC), yc, xc, ops.whiteningRange);
else
    %
    [E, D] 	= svd(CC);
    D = diag(D);
    eps 	= 1e-6;
    Wrot 	= E * diag(1./(D + eps).^.5) * E';
end
Wrot    = ops.scaleproc * Wrot;

fprintf('Time %3.0fs. Loading raw data and applying filters... \n', toc);

fid         = fopen(ops.fbinary, 'r');
fidW    = fopen(ops.fproc, 'w');

if strcmp(ops.initialize, 'fromData')
    i0  = 0;
    ixt  = round(linspace(1, size(ops.wPCA,1), ops.nt0));
    wPCA = ops.wPCA(ixt, 1:3);
    
    rez.ops.wPCA = wPCA; % write wPCA back into the rez structure
    uproj = zeros(1e6,  size(wPCA,2) * Nchan, 'single');
end
%
for ibatch = 1:Nbatch
    if isproc(ibatch) %ibatch<=Nbatch_buff
        if ops.GPU
            datr = single(gpuArray(DATA(:,:,ibatch)));
        else
            datr = single(DATA(:,:,ibatch));
        end
    else
        offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
        if ibatch==1
            ioffset = 0;
        else
            ioffset = ops.ntbuff;
        end
        fseek(fid, offset, 'bof');
        
        buff = fread(fid, [NchanTOT NTbuff], '*int16');
        if isempty(buff)
            break;
        end
        nsampcurr = size(buff,2);
        if nsampcurr<NTbuff
            buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
        end
        
        if ops.GPU
            dataRAW = gpuArray(buff);
        else
            dataRAW = buff;
        end
        dataRAW = dataRAW';
        dataRAW = single(dataRAW);
        dataRAW = dataRAW(:, chanMapConn);
        
        datr = filter(b1, a1, dataRAW);
        datr = flipud(datr);
        datr = filter(b1, a1, datr);
        datr = flipud(datr);
        
        datr = datr(ioffset + (1:NT),:);
    end
    
    datr    = datr * Wrot;
    
    if ops.GPU
        dataRAW = gpuArray(datr);
    else
        dataRAW = datr;
    end
    %         dataRAW = datr;
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    if strcmp(ops.initialize, 'fromData') %&& rem(ibatch, 10)==1
        % find isolated spikes
        [row, col, mu] = isolated_peaks_(dataRAW, ops.loc_range, ops.long_range, ops.spkTh);
        
        % find their PC projections
        uS = get_PCproj_(dataRAW, row, col, wPCA, ops.maskMaxChannels);
        
        uS = permute(uS, [2 1 3]);
        uS = reshape(uS,numel(row), Nchan * size(wPCA,2));
        
        if i0+numel(row)>size(uproj,1)
            uproj(1e6 + size(uproj,1), 1) = 0;
        end
        
        uproj(i0 + (1:numel(row)), :) = gather_try_(uS);
        i0 = i0 + numel(row);
    end
    
    if ibatch<=Nbatch_buff
        DATA(:,:,ibatch) = gather_try_(datr);
    else
        datcpu  = gather_try_(int16(datr));
        fwrite(fidW, datcpu, 'int16');
    end
    
end

if strcmp(ops.initialize, 'fromData')
   uproj(i0+1:end, :) = []; 
end
Wrot        = gather_try_(Wrot);
rez.Wrot    = Wrot;

fclose(fidW);
fclose(fid);
if ops.verbose
    fprintf('Time %3.2f. Whitened data written to disk... \n', toc);
    fprintf('Time %3.2f. Preprocessing complete!\n', toc);
end


rez.temp.Nbatch = Nbatch;
rez.temp.Nbatch_buff = Nbatch_buff;
end %func


%--------------------------------------------------------------------------
function rez = fitTemplates_(rez, DATA, uproj)
nt0             = rez.ops.nt0;
rez.ops.nt0min  = ceil(20 * nt0/61);

ops = rez.ops;

rng('default');
rng(1);

Nbatch      = rez.temp.Nbatch;
Nbatch_buff = rez.temp.Nbatch_buff;

Nfilt 	= ops.Nfilt; %256+128;

ntbuff  = ops.ntbuff;
NT  	= ops.NT;

Nrank   = ops.Nrank;
Th 		= ops.Th;
maxFR 	= ops.maxFR;

Nchan 	= ops.Nchan;

batchstart = 0:NT:NT*(Nbatch-Nbatch_buff);

delta = NaN * ones(Nbatch, 1);
iperm = randperm(Nbatch);

switch ops.initialize
    case 'fromData'
        WUinit = optimizePeaks_(ops,uproj);%does a scaled kmeans 
        dWU    = WUinit(:,:,1:Nfilt);
        %             dWU = alignWU(dWU);
    otherwise
        initialize_waves0;
        ipck = randperm(size(Winit,2), Nfilt);
        W = [];
        U = [];
        for i = 1:Nrank
            W = cat(3, W, Winit(:, ipck)/Nrank);
            U = cat(3, U, Uinit(:, ipck));
        end
        W = alignW_(W, ops);
        
        dWU = zeros(nt0, Nchan, Nfilt, 'single');
        for k = 1:Nfilt
            wu = squeeze(W(:,k,:)) * squeeze(U(:,k,:))';
            newnorm = sum(wu(:).^2).^.5;
            W(:,k,:) = W(:,k,:)/newnorm;
            
            dWU(:,:,k) = 10 * wu;
        end
        WUinit = dWU;
end
[W, U, mu, UtU, nu] = decompose_dWU_(ops, dWU, Nrank, rez.ops.kcoords);
W0 = W;
W0(NT, 1) = 0;
fW = fft(W0, [], 1);
fW = conj(fW);

nspikes = zeros(Nfilt, Nbatch);
lam =  ones(Nfilt, 1, 'single');

freqUpdate = 100 * 4;
iUpdate = 1:freqUpdate:Nbatch;


dbins = zeros(100, Nfilt);
dsum = 0;
miniorder = repmat(iperm, 1, ops.nfullpasses);
%     miniorder = repmat([1:Nbatch Nbatch:-1:1], 1, ops.nfullpasses/2);

i = 1; % first iteration

epu = ops.epu;


%-----
% pmi = exp(-1./exp(linspace(log(ops.momentum(1)), log(ops.momentum(2)), Nbatch*ops.nannealpasses)));
pmi = exp(-1./linspace(1/ops.momentum(1), 1/ops.momentum(2), Nbatch*ops.nannealpasses));
% pmi = exp(-linspace(ops.momentum(1), ops.momentum(2), Nbatch*ops.nannealpasses));

% pmi  = linspace(ops.momentum(1), ops.momentum(2), Nbatch*ops.nannealpasses);
Thi  = linspace(ops.Th(1),                 ops.Th(2), Nbatch*ops.nannealpasses);
if ops.lam(1)==0
    lami = linspace(ops.lam(1), ops.lam(2), Nbatch*ops.nannealpasses);
else
    lami = exp(linspace(log(ops.lam(1)), log(ops.lam(2)), Nbatch*ops.nannealpasses));
end

if Nbatch_buff<Nbatch
    fid = fopen(ops.fproc, 'r');
end

nswitch = [0];
msg = [];
fprintf('Time %3.0fs. Optimizing templates ...\n', toc)
while (i<=Nbatch * ops.nfullpasses+1)
    % set the annealing parameters
    if i<Nbatch*ops.nannealpasses
        Th      = Thi(i);
        lam(:)  = lami(i);
        pm      = pmi(i);
    end
    
    % some of the parameters change with iteration number
    Params = double([NT Nfilt Th maxFR 10 Nchan Nrank pm epu nt0]);
    
    % update the parameters every freqUpdate iterations
    if i>1 &&  ismember(rem(i,Nbatch), iUpdate) %&& i>Nbatch
        dWU = gather_try_(dWU);
        
        % break bimodal clusters and remove low variance clusters
        if  ops.shuffle_clusters &&...
                i>Nbatch && rem(rem(i,Nbatch), 4*400)==1    % i<Nbatch*ops.nannealpasses
            [dWU, dbins, nswitch, nspikes, iswitch] = ...
                replace_clusters_(dWU, dbins,  Nbatch, ops.mergeT, ops.splitT, WUinit, nspikes);
        end
        
        dWU = alignWU_(dWU, ops);
        
        % restrict spikes to their peak group
        %         dWU = decompose_dWU(dWU, kcoords);
        
        % parameter update
        [W, U, mu, UtU, nu] = decompose_dWU_(ops, dWU, Nrank, rez.ops.kcoords);
        
        if ops.GPU
            dWU = gpuArray(dWU);
        else
            W0 = W;
            W0(NT, 1) = 0;
            fW = fft(W0, [], 1);
            fW = conj(fW);
        end
        
        NSP = sum(nspikes,2);
        if ops.showfigures
%             set(0,'DefaultFigureWindowStyle','docked')
%             figure;
            subplot(2,2,1)
            for j = 1:10:Nfilt
                if j+9>Nfilt;
                    j = Nfilt -9;
                end
                plot(log(1+NSP(j + [0:1:9])), mu(j+ [0:1:9]), 'o');
                xlabel('log of number of spikes')
                ylabel('amplitude of template')
                hold all
            end
            axis tight;
            title(sprintf('%d  ', nswitch));
            subplot(2,2,2)
            plot(W(:,:,1))
            title('timecourses of top PC')
            
            subplot(2,2,3)
            imagesc(U(:,:,1))
            title('spatial mask of top PC')
            
            drawnow
        end
        % break if last iteration reached
        if i>Nbatch * ops.nfullpasses; break; end
        
        % record the error function for this iteration
        rez.errall(ceil(i/freqUpdate))          = nanmean(delta);
        
    end
    
    % select batch and load from RAM or disk
    ibatch = miniorder(i);
    if ibatch>Nbatch_buff
        offset = 2 * ops.Nchan*batchstart(ibatch-Nbatch_buff);
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NT ops.Nchan], '*int16');
    else
        dat = DATA(:,:,ibatch);
    end
    
    % move data to GPU and scale it
    if ops.GPU
        dataRAW = gpuArray(dat);
    else
        dataRAW = dat;
    end
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    % project data in low-dim space
    data = dataRAW * U(:,:);
    
    if ops.GPU
        % run GPU code to get spike times and coefficients
        [dWU, ~, id, x,Cost, nsp] = ...
            mexMPregMU(Params,dataRAW,W,data,UtU,mu, lam .* (20./mu).^2, dWU, nu);
    else
        [dWU, ~, id, x,Cost, nsp] = ...
            mexMPregMUcpu(Params,dataRAW,fW,data,UtU,mu, lam .* (20./mu).^2, dWU, nu, ops);
    end
    
    dbins = .9975 * dbins;  % this is a hard-coded forgetting factor, needs to become an option
    if ~isempty(id)
        % compute numbers of spikes
        nsp                = gather_try_(nsp(:));
        nspikes(:, ibatch) = nsp;
        
        % bin the amplitudes of the spikes
        xround = min(max(1, int32(x)), 100);
        
        dbins(xround + id * size(dbins,1)) = dbins(xround + id * size(dbins,1)) + 1;
        
        % estimate cost function at this time step
        delta(ibatch) = sum(Cost)/1e3;
    end
    
    % update status
    if ops.verbose  && rem(i,20)==1
        nsort = sort(round(sum(nspikes,2)), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Time %2.2f, batch %d/%d, mu %2.2f, neg-err %2.6f, NTOT %d, n100 %d, n200 %d, n300 %d, n400 %d\n', ...
            toc, i,Nbatch* ops.nfullpasses,nanmean(mu(:)), nanmean(delta), round(sum(nsort)), ...
            nsort(min(size(W,2), 100)), nsort(min(size(W,2), 200)), ...
            nsort(min(size(W,2), 300)), nsort(min(size(W,2), 400)));
        fprintf(msg);
    end
    
    % increase iteration counter
    i = i+1;
end

% close the data file if it has been used
if Nbatch_buff<Nbatch
    fclose(fid);
end

if ~ops.GPU
   rez.fW = fW; % save fourier space templates if on CPU
end
rez.dWU               = gather_try_(dWU);
rez.nspikes               = nspikes;
end %func


%--------------------------------------------------------------------------
function rez = fullMPMU_(rez, DATA)
ops = rez.ops;

Nfilt   = ops.Nfilt;
lam =  ones(Nfilt, 1, 'single');
lam(:)    = ops.lam(3);

[W, U, mu, UtU, nu] = decompose_dWU_(ops, rez.dWU, ops.Nrank,rez.ops.kcoords);


pm = exp(-ops.momentum(2));
Params = double([ops.NT ops.Nfilt ops.Th(3) ops.maxFR 10 ops.Nchan ops.Nrank pm ops.epu ops.nt0]);

Params(3) = ops.Th(3);
Params(4) = 50000;
Params(5) = 50; 

if ops.GPU
    U0 = gpuArray(U);
else
    U0 = U;
end


%-----
nt0     = rez.ops.nt0;
Nrank   = ops.Nrank;
WtW     = zeros(Nfilt,Nfilt,2*nt0-1, 'single');
for i = 1:Nrank
    for j = 1:Nrank
        utu0 = U0(:,:,i)' * U0(:,:,j);
        if ops.GPU
            wtw0 =  gather_try_(mexWtW2(Params, W(:,:,i), W(:,:,j), utu0));
        else
            wtw0 =  getWtW2_(Params, W(:,:,i), W(:,:,j), utu0);
            wtw0 = permute(wtw0, [2 3 1]);
        end
        WtW = WtW + wtw0;
        clear wtw0 utu0
%         wtw0 = squeeze(wtw(:,i,:,j,:));
        
    end
end

mWtW = max(WtW, [], 3);
WtW = permute(WtW, [3 1 2]);

if ops.GPU
    WtW = gpuArray(WtW);
end


%-----
Nbatch_buff = rez.temp.Nbatch_buff;
Nbatch      = rez.temp.Nbatch;
Nchan       = ops.Nchan;
if ~ops.GPU
   fW = rez.fW; % load fft-ed templates 
end
% mWtW = mWtW - diag(diag(mWtW));

% rez.WtW = gather_try_(WtW);
%
clear wtw0 utu0 U0
%
clear nspikes2
st3 = [];
rez.st3 = [];

if ops.verbose
   fprintf('Time %3.0fs. Running the final template matching pass...\n', toc) 
end

if Nbatch_buff<Nbatch
    fid = fopen(ops.fproc, 'r');
end
msg = [];

if ~isempty(ops.nNeigh)
    nNeigh    = ops.nNeigh;
    
    rez.cProj = zeros(5e6, nNeigh, 'single');

    % sort pairwise templates
    nsp = sum(rez.nspikes,2);
    vld = single(nsp>100);
    cr    = mWtW .* (vld * vld');
    cr(isnan(cr)) = 0;
    [~, iNgsort] = sort(cr, 1, 'descend');
    
    % save full similarity score
    rez.simScore = cr;
    maskTT = zeros(Nfilt, 'single');
    rez.iNeigh = iNgsort(1:nNeigh, :);
    for i = 1:Nfilt
        maskTT(rez.iNeigh(:,i),i) = 1;
    end
end
if ~isempty(ops.nNeighPC)
    nNeighPC  = ops.nNeighPC;
    load('./kilosort/PCspikes.mat');
    ixt = round(linspace(1, size(Wi,1), ops.nt0));
    Wi = Wi(ixt, 1:3);
    rez.cProjPC = zeros(5e6, 3*nNeighPC, 'single');
    
    % sort best channels
    [~, iNch]       = sort(abs(U(:,:,1)), 1, 'descend');
    maskPC          = zeros(Nchan, Nfilt, 'single');
    rez.iNeighPC    = iNch(1:nNeighPC, :);
    for i = 1:Nfilt
        maskPC(rez.iNeighPC(:,i),i) = 1;
    end
    maskPC = repmat(maskPC, 3, 1);
end

irun = 0;
i1nt0 = int32([1:nt0])';


%-----
LAM = lam .* (20./mu).^2;

NT = ops.NT;
batchstart = 0:NT:NT*(Nbatch-Nbatch_buff);

for ibatch = 1:Nbatch    
    if ibatch>Nbatch_buff
        offset = 2 * ops.Nchan*batchstart(ibatch-Nbatch_buff); % - ioffset;
        fseek(fid, offset, 'bof');
        dat = fread(fid, [NT ops.Nchan], '*int16');
    else
       dat = DATA(:,:,ibatch); 
    end
    if ops.GPU
        dataRAW = gpuArray(dat);
    else
        dataRAW = dat;
    end
    dataRAW = single(dataRAW);
    dataRAW = dataRAW / ops.scaleproc;
    
    % project data in low-dim space
    if ops.GPU
        data    = gpuArray.zeros(NT, Nfilt, Nrank, 'single');
    else
        data   = zeros(NT, Nfilt, Nrank, 'single');
    end
    for irank = 1:Nrank
        data(:,:,irank) = dataRAW * U(:,:,irank);
    end
    data                = reshape(data, NT, Nfilt*Nrank);

    if ops.GPU
        [st, id, x, errC, PCproj] ...
                        = mexMPmuFEAT(Params,data,W,WtW, mu, lam .* (20./mu).^2, nu);
    else
         [st, id, x, errC, PCproj]= cpuMPmuFEAT_(Params,data,fW,WtW, mu, lam .* (20./mu).^2, nu, ops);
    end
    
    if ~isempty(st)
        if ~isempty(ops.nNeighPC)
            % PCA coefficients
            inds            = repmat(st', nt0, 1) + repmat(i1nt0, 1, numel(st));
            try  datSp      = dataRAW(inds(:), :);
            catch
                datSp       = dataRAW(inds(:), :);
            end
            datSp           = reshape(datSp, [size(inds) Nchan]);
            coefs           = reshape(Wi' * reshape(datSp, nt0, []), size(Wi,2), numel(st), Nchan);
            coefs           = reshape(permute(coefs, [3 1 2]), [], numel(st));
            coefs           = coefs .* maskPC(:, id+1);
            iCoefs          = reshape(find(maskPC(:, id+1)>0), 3*nNeighPC, []);
            rez.cProjPC(irun + (1:numel(st)), :) = gather_try_(coefs(iCoefs)');
        end
        if ~isempty(ops.nNeigh)
            % template coefficients
            % transform coefficients
            PCproj          = bsxfun(@rdivide, ...
                bsxfun(@plus, PCproj, LAM.*mu), sqrt(1+LAM));
            
            PCproj          = maskTT(:, id+1) .* PCproj;
            iPP             = reshape(find(maskTT(:, id+1)>0), nNeigh, []);
            rez.cProj(irun + (1:numel(st)), :) = PCproj(iPP)';
        end
        % increment number of spikes
        irun            = irun + numel(st);
        
        if ibatch==1;
            ioffset         = 0;
        else
            ioffset         = ops.ntbuff;
        end
        st                  = st - ioffset;
        
        %     nspikes2(1:size(W,2)+1, ibatch) = histc(id, 0:1:size(W,2));
        STT = cat(2, ops.nt0min + double(st) +(NT-ops.ntbuff)*(ibatch-1), ...
            double(id)+1, double(x), ibatch*ones(numel(x),1));
        st3             = cat(1, st3, STT);
    end
    if rem(ibatch,100)==1
%         nsort = sort(sum(nspikes2,2), 'descend');
        fprintf(repmat('\b', 1, numel(msg)));
        msg             = sprintf('Time %2.2f, batch %d/%d,  NTOT %d\n', ...
            toc, ibatch,Nbatch, size(st3,1));        
        fprintf(msg);        
    end
end


%-----
[~, isort]      = sort(st3(:,1), 'ascend');
st3             = st3(isort,:);

rez.st3         = st3;
if ~isempty(ops.nNeighPC)
    % re-sort coefficients for projections
    rez.cProjPC(irun+1:end, :)  = [];
    rez.cProjPC                 = reshape(rez.cProjPC, size(rez.cProjPC,1), [], 3);
    rez.cProjPC                 = rez.cProjPC(isort, :,:);
    for ik = 1:Nfilt
        iSp                     = rez.st3(:,2)==ik;
        OneToN                  = 1:nNeighPC;
        [~, isortNeigh]         = sort(rez.iNeighPC(:,ik), 'ascend');
        OneToN(isortNeigh)      = OneToN;
        rez.cProjPC(iSp, :,:)   = rez.cProjPC(iSp, OneToN, :);
    end
    
    rez.cProjPC                 = permute(rez.cProjPC, [1 3 2]);
end
if ~isempty(ops.nNeigh)
    rez.cProj(irun+1:end, :)    = [];
    rez.cProj                   = rez.cProj(isort, :);

    % re-index the template coefficients
    for ik = 1:Nfilt
        iSp                     = rez.st3(:,2)==ik;
        OneToN                  = 1:nNeigh;
        [~, isortNeigh]         = sort(rez.iNeigh(:,ik), 'ascend');
        OneToN(isortNeigh)      = OneToN;
        rez.cProj(iSp, :)       = rez.cProj(iSp, OneToN);
    end
end


%-----
% rez.ops             = ops;
rez.W               = W;
rez.U               = U;
rez.mu              = mu;

rez.t2p = [];
for i = 1:Nfilt
    wav0            = W(:,i,1);
    wav0            = my_conv_(wav0', .5)';
   [~, itrough]     = min(wav0);
    [~, t2p]        = max(wav0(itrough:end));
    rez.t2p(i,1)    = t2p;
    rez.t2p(i,2)    = itrough;   
end

rez.nbins           = histc(rez.st3(:,2), .5:1:Nfilt+1);

[~, rez.ypos]       = max(rez.U(:,:,1), [], 1);
if Nbatch_buff<Nbatch
    fclose(fid);
end

% center the templates
rez.W               = cat(1, zeros(nt0 - (ops.nt0-1-ops.nt0min), Nfilt, Nrank), rez.W);
rez.WrotInv         = (rez.Wrot/200)^-1;


%-----
Urot = U;
for k = 1:size(U,3)
   Urot(:,:,k)  = rez.WrotInv' * Urot(:,:,k);
end
for n = 1:size(U,2)
    rez.Wraw(:,:,n) = mu(n) * sq_(Urot(:,n,:)) * sq_(rez.W(:,n,:))';
end
end %func


%--------------------------------------------------------------------------
function ops = convertOpenEphysToRawBInary_(ops)

fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary)); 
fidout      = fopen(fname, 'w');
%
clear fs
for j = 1:ops.Nchan
   fs{j} = dir(fullfile(ops.root, sprintf('*CH%d_*.continuous', j) ));
end
nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
%
nBlocks     = unique(nblocks);
nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, 1);

tic
for k = 1:nBlocks
    for j = 1:ops.Nchan
        fid{j}             = fopen(fullfile(ops.root, fs{j}(k).name));
        % discard header information
        fseek(fid{j}, 1024, 0);
    end
    %
    nsamps = 0;
    flag = 1;
    while 1
        samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
        for j = 1:ops.Nchan
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');
            
            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');

            nbatches        = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end
        
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
       
        samples         = samples';
        fwrite(fidout, samples, 'int16');
        
        nsamps = nsamps + size(samples,2);
        
        if flag==0
            break;
        end
    end
    ops.nSamplesBlocks(k) = nsamps;
    
    for j = 1:ops.Nchan
       fclose(fid{j}); 
    end
    
end
    
fclose(fidout);

toc
end %func


%--------------------------------------------------------------------------
function x = gather_try_(x)
try
    x = gather(x);
catch
end
end %func


%--------------------------------------------------------------------------
function v = getOr_(s, field, default)
%getOr Returns the structure field or a default if either don't exist
%   v = getOr(s, field, [default]) returns the 'field' of the structure 's'
%   or 'default' if the structure is empty of the field does not exist. If
%   default is not specified it defaults to []. 'field' can also be a cell
%   array of fields, in which case it will check for all of them and return
%   the value of the first field that exists, if any (otherwise the default
%   value).

if nargin < 3
  default = [];
end

fieldExists = isfield(s, field);
if any(fieldExists)
  if iscellstr(field)
    v = s.(field{find(fieldExists, 1)});
  else
    v = s.(field);
  end
else
  v = default;
end
end %func


%--------------------------------------------------------------------------
function Us = get_PCproj_(S1, row, col, wPCA, maskMaxChans)

[nT, nChan] = size(S1);
dt = -21 + [1:size(wPCA,1)];
inds = repmat(row', numel(dt), 1) + repmat(dt', 1, numel(row));

clips = reshape(S1(inds, :), numel(dt), numel(row), nChan);


mask = repmat([1:nChan], [numel(row) 1]) - repmat(col, 1, nChan);
Mask(1,:,:) = abs(mask)<maskMaxChans;

clips = bsxfun(@times, clips , Mask);

Us = wPCA' * reshape(clips, numel(dt), []);
Us = reshape(Us, size(wPCA,2), numel(row), nChan);

Us = permute(Us, [3 2 1]);
end %func


%--------------------------------------------------------------------------
function [row, col, mu] = isolated_peaks_(S1, loc_range, long_range, Th)
% loc_range = [3  1];
% long_range = [30  6];
smin = my_min_(S1, loc_range, [1 2]);
peaks = single(S1<smin+1e-3 & S1<Th);

sum_peaks = my_sum_(peaks, long_range, [1 2]);
peaks = peaks .* (sum_peaks<1.2).* S1;

peaks([1:20 end-40:end], :) = 0;

[row, col, mu] = find(peaks);
end %func


%--------------------------------------------------------------------------
function S1 = my_min_(S1, sig, varargin)
% takes an extra argument which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end
if numel(idims)>1 && numel(sig)>1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

for i = 1:length(idims)
    sig = sigall(i);
    
    idim = idims(i);
    Nd = ndims(S1);
    
    S1 = permute(S1, [idim 1:idim-1 idim+1:Nd]);
    
    dsnew = size(S1);

    S1 = reshape(S1, size(S1,1), []);
    dsnew2 = size(S1);
    
    S1 = cat(1, Inf*ones([sig, dsnew2(2)]),S1, Inf*ones([sig, dsnew2(2)]));
    Smax = S1(1:dsnew2(1), :);
    for j = 1:2*sig
        Smax = min(Smax, S1(j + (1:dsnew2(1)), :));
    end
    
    S1 = reshape(Smax, dsnew);
       
    S1 = permute(S1, [2:idim 1 idim+1:Nd]);
end
end %func


%--------------------------------------------------------------------------
function Wrot = whiteningLocal_(CC, yc, xc, nRange)

Wrot = zeros(size(CC,1), size(CC,1));
for j = 1:size(CC,1)
    ds          = (xc - xc(j)).^2 + (yc - yc(j)).^2;
    [~, ilocal] = sort(ds, 'ascend');
    ilocal      = ilocal(1:nRange);
    
    [E, D]      = svd(CC(ilocal, ilocal));
    D           = diag(D);
    eps         = 1e-6;
    wrot0       = E * diag(1./(D + eps).^.5) * E';
    Wrot(ilocal, j)  = wrot0(:,1);
end
end %func


%--------------------------------------------------------------------------
function W = alignW_(W, ops)

[nt0 , Nfilt] = size(W);


[~, imax] = min(W, [], 1);
dmax = -(imax - ops.nt0min);
% dmax = min(1, abs(dmax)) .* sign(dmax);
 
for i = 1:Nfilt
    if dmax(i)>0
        W((dmax(i) + 1):nt0, i) = W(1:nt0-dmax(i), i);
    else
        W(1:nt0+dmax(i), i) = W((1-dmax(i)):nt0, i);
    end
end
end %func


%--------------------------------------------------------------------------
function WU = alignWU_(WU, ops)

[nt0 , Nchan, Nfilt] = size(WU);
[~, imin] = min(reshape(WU, nt0*Nchan, Nfilt), [], 1);

iMinChan = ceil(imin/nt0);


% imin = rem(imin-1, nt0) + 1;

% [~, imax] = min(W, [], 1);
% dmax = -(imin - 20);
% dmax = min(1, abs(dmax)) .* sign(dmax);
 
dmax = zeros(Nfilt, 1);
for i = 1:Nfilt
    wu = WU(:,iMinChan(i),i);
%     [~, imin] = min(diff(wu, 1));
    [~, imin] = min(wu);
    dmax(i) = - (imin- ops.nt0min);
    
    if dmax(i)>0
        WU((dmax(i) + 1):nt0, :,i) = WU(1:nt0-dmax(i),:, i);
    else
        WU(1:nt0+dmax(i),:, i) = WU((1-dmax(i)):nt0,:, i);
    end
end
end %func


%--------------------------------------------------------------------------
function [W, U, mu, UtU, nu] = decompose_dWU_(ops, dWU, Nrank, kcoords)

[nt0 Nchan Nfilt] = size(dWU);

W = zeros(nt0, Nrank, Nfilt, 'single');
U = zeros(Nchan, Nrank, Nfilt, 'single');
mu = zeros(Nfilt, 1, 'single');
% dmax = zeros(Nfilt, 1);

dWU(isnan(dWU)) = 0;
if ops.parfor
    parfor k = 1:Nfilt
        [W(:,:,k), U(:,:,k), mu(k)] = get_svds_(dWU(:,:,k), Nrank);
    end
else
    for k = 1:Nfilt
        [W(:,:,k), U(:,:,k), mu(k)] = get_svds_(dWU(:,:,k), Nrank);
    end
end
U = permute(U, [1 3 2]);
W = permute(W, [1 3 2]);

U(isnan(U)) = 0;

if numel(unique(kcoords))>1
    U = zeroOutKcoords_(U, kcoords, ops.criterionNoiseChannels);
end

UtU = abs(U(:,:,1)' * U(:,:,1)) > .1;


Wdiff = cat(1, W, zeros(2, Nfilt, Nrank)) - cat(1,  zeros(2, Nfilt, Nrank), W);
nu = sum(sum(Wdiff.^2,1),3);
nu = nu(:);



% mu = min(mu, 200);
end %func


%--------------------------------------------------------------------------
% addpath('C:\CODE\GitHub\KiloSort\preDetect')
function WUinit = optimizePeaks_(ops,uproj)

nt0             = ops.nt0;

nProj           = size(uproj,2);
nSpikesPerBatch = 4000;
inds            = 1:nSpikesPerBatch * floor(size(uproj,1)/nSpikesPerBatch);
inds            = reshape(inds, nSpikesPerBatch, []);
% Nbatch = size(inds,2);
iperm           = randperm(size(inds,2));
miniorder       = repmat(iperm, 1, ops.nfullpasses);
%     miniorder = repmat([1:Nbatch Nbatch:-1:1], 1, ops.nfullpasses/2);

if ~exist('spikes_merged', 'var')
    uBase = zeros(1e4, nProj);
    nS = zeros(1e4, 1);
    ncurr = 1;
    
    for ibatch = 1:size(inds,2)
        % merge in with existing templates
        uS = uproj(inds(:,ibatch), :);
        [nSnew, iNonMatch] = merge_spikes0_(uBase(1:ncurr,:), nS(1:ncurr), uS, ops.crit);
        nS(1:ncurr) = nSnew;
        %
        % reduce non-matches
        [uNew, nSadd] = reduce_clusters0_(uS(iNonMatch,:), ops.crit);
        
        % add new spikes to list
        uBase(ncurr + [1:size(uNew,1)], :) = uNew;
        nS(ncurr + [1:size(uNew,1)]) = nSadd;
        
        ncurr = ncurr + size(uNew,1);
        
        if ncurr>1e4
            break;
        end
    end
    %
    nS = nS(1:ncurr);
    uBase = uBase(1:ncurr, :);
    spikes_merged = 1;
end
[~, itsort] = sort(nS, 'descend');

%-----
% initialize U
Nfilt = ops.Nfilt;
lam = ops.lam(1) * ones(Nfilt, 1, 'single');

ind_filt = itsort(rem([1:Nfilt]-1, numel(itsort)) + 1);
U = gpuArray(uBase(ind_filt, :))';
U = U + .001 * randn(size(U));
mu = sum(U.^2,1)'.^.5;
U = normc_(U);

for i = 1:10
    
    idT = zeros(size(inds));
    dWU = zeros(Nfilt, nProj, 'single');
    nToT = gpuArray.zeros(Nfilt, 1, 'single');
    Cost = gpuArray(single(0));
    
    for ibatch = 1:size(inds,2)
        % find clusters
        clips = reshape(gpuArray(uproj(inds(:,ibatch), :)), nSpikesPerBatch, nProj);
        
        ci = clips * U;
        
        ci = bsxfun(@plus, ci, (mu .* lam)');
        cf = bsxfun(@rdivide, ci.^2, 1 + lam');
        cf = bsxfun(@minus, cf, (mu.^2.*lam)');
        
        [max_cf, id] = max(cf, [], 2);
        
        id = gather_try_(id);
        %        x = ci([1:nSpikesPerBatch] + nSpikesPerBatch * (id-1)')' - mu(id) .* lam(id);
        idT(:,ibatch) = id;
        
        L = gpuArray.zeros(Nfilt, nSpikesPerBatch, 'single');
        L(id' + [0:Nfilt:(Nfilt*nSpikesPerBatch-1)]) = 1;
        dWU = dWU + L * clips;
        
        nToT = nToT + sum(L, 2);
        Cost = Cost + mean(max_cf);
    end
    dWU  = bsxfun(@rdivide, dWU, nToT);
    
    U = dWU';
    mu = sum(U.^2,1)'.^.5;
    U = normc_(U);
    Cost = Cost/size(inds,2);
    
%     disp(Cost)
    
%     plot(sort(log(1+nToT)))
%     drawnow
end

%-----
Nchan = ops.Nchan;
Nfilt = ops.Nfilt;
wPCA = ops.wPCA(:,1:3);
Urec = reshape(U, Nchan, size(wPCA,2), Nfilt);

Urec= permute(Urec, [2 1 3]);
Wrec = reshape(wPCA * Urec(:,:), nt0, Nchan, Nfilt);

Wrec = gather_try_(Wrec);
Nrank = 3;
W = zeros(nt0, Nfilt, Nrank, 'single');
U = zeros(Nchan, Nfilt, Nrank, 'single');

Wrec(isnan(Wrec(:))) = 0;
for j = 1:Nfilt
    [w sv u] = svd(Wrec(:,:,j));
    w = w * sv;
    
    Sv = diag(sv);
    W(:,j,:) = w(:, 1:Nrank)/sum(Sv(1:ops.Nrank).^2).^.5;
    U(:,j,:) = u(:, 1:Nrank);
end

Uinit = U;
Winit = W;
mu = gather_try_(single(mu));
muinit = mu;

WUinit = zeros(nt0, Nchan, Nfilt);
for j = 1:Nfilt
    WUinit(:,:,j) = muinit(j)  * Wrec(:,:,j);
end
WUinit = single(WUinit);
end %func


%--------------------------------------------------------------------------
function [dWUtot, dbins, nswitch, nspikes, iYout] = ...
    replace_clusters_(dWUtot,dbins, Nbatch, mergeT, splitT,  WUinit, nspikes)

uu = Nbatch * dbins;
nhist = 1:1:100;
% nSpikes = sum(uu,1);
nSpikes = sum(nspikes,2)';

[score, iY1, mu1, mu2, u1, u2]   = split_clust_(uu, nhist);
[~, iY, drez]                    = distance_betwxt_(dWUtot);

[dsort, isort] = sort(drez, 'ascend');
iYout = iY(isort);

nmerged = sum(dsort<mergeT);
nsplit = sum(score>splitT);

mu = sum(sum(dWUtot.^2,1),2).^.5;
mu = mu(:);
freeInd = find(nSpikes<200 | mu'<10 | isnan(mu'));

for k = 1:nmerged
    % merge the two clusters
    iMerged = iY(isort(k));
    wt = [nSpikes(iMerged); nSpikes(isort(k))];
    wt = wt/sum(wt);
%     mu(iMerged) = [mu(iMerged) mu(isort(k))] * wt;
    
    dWUtot(:,:,iMerged)  = dWUtot(:,:,iMerged) * wt(1) + dWUtot(:,:,isort(k)) * wt(2);
    dWUtot(:,:,isort(k)) = 1e-10;
    
    nspikes(iMerged, :) = nspikes(iMerged, :) + nspikes(isort(k), :);
    nspikes(isort(k), :) = 0;
end


for k = 1:min(nmerged+numel(freeInd), nsplit)
    if k<=numel(freeInd)
        inew= freeInd(k);
    else
        inew = isort(k - numel(freeInd));
    end
    
    mu0 = mu(iY1(k));
    
    % split the bimodal cluster, overwrite merged cluster
    mu(inew)     = mu1(k);
    mu(iY1(k))   = mu2(k);
    
    dbins(:, inew)     = u1(:, k) /Nbatch;
    dbins(:, iY1(k))   = u2(:, k) /Nbatch;

    nspikes(inew, :)     = nspikes(iY1(k), :)/2;
    nspikes(iY1(k), :)   = nspikes(iY1(k), :)/2;
    dWUtot(:,:,inew)     = mu1(k)/mu0 * dWUtot(:,:,iY1(k)); %/npm(iY1(k));
    dWUtot(:,:,iY1(k))   = mu2(k)/mu0 * dWUtot(:,:,iY1(k)); %/npm(iY1(k));
end

d2d                 = pairwise_dists_(dWUtot, WUinit);
dmatch              = min(d2d, [], 1);

[~, inovel] = sort(dmatch, 'descend');
% inovel = find(dmatch(1:1000)>.4);
% inovel = inovel(randperm(numel(inovel)));

i0 = 0;

for k = 1+min(nmerged+numel(freeInd), nsplit):nmerged+numel(freeInd)
    % add new clusters
    i0 = i0 + 1;
    if i0>numel(inovel)
        break;
    end
    if k<=numel(freeInd)
        inew= freeInd(k);
    else
        inew = isort(k - numel(freeInd));
    end
     
    dbins(:, inew)     = 1;
    
    nspikes(inew, :) = 1/8;
    
    
    dWUtot(:,:,inew)     = WUinit(:,:,inovel(i0)); %ratio * mu1(k)/mu0 * dWUtot(:,:,iY1(k));
    
end

nswitch = [min(nmerged, nsplit) i0]; %min(nmerged+numel(freeInd), nsplit);
end %func


%--------------------------------------------------------------------------
function [sts, ids, xs, Costs, cprojall] = cpuMPmuFEAT_(Params,data,fW,WtW, mu, lam1, nu, ops)

nt0 = ops.nt0;

WtW     = permute(WtW, [1 3 2]);

NT      = Params(1);
nFilt   = Params(2);
Th      = Params(3);

fdata   = fft(data, [], 1);
proj    = real(ifft(fdata .* fW(:,:), [], 1));
proj    = sum(reshape(proj, NT, nFilt, 3),3);

trange = int32([-(nt0-1):(nt0-1)]);

xs      = zeros(Params(4), 1, 'single');
ids     = zeros(Params(4), 1, 'int32');
sts     = zeros(Params(4), 1, 'int32');
Costs    = zeros(Params(4), 1, 'single');
cprojall = zeros(Params(4), nFilt, 'single');

i0 = 0;
for k = 1:30
    Ci = bsxfun(@plus, proj, (mu.*lam1)');
    Ci = bsxfun(@rdivide, Ci.^2,  1 + lam1');
    Ci = bsxfun(@minus, Ci, (lam1 .* mu.^2)');
    
    [mX, id] = max(Ci,[], 2);
    
    maX         = -my_min_(-mX, 31, 1);
    id          = int32(id);
    
    st                   = find((maX < mX + 1e-3) & mX > Th*Th);
    st(st>NT-nt0 | st<nt0) = [];
    
    if isempty(st)
       break; 
    end
    id      = id(st);
    
    % inds = bsxfun(@plus, st', [1:nt0]');
    
    x       = zeros(size(id));
    Cost    = zeros(size(id));
    nsp     = zeros(nFilt,1);
    cproj   = zeros(size(id,1), nFilt, 'single');
    for j = 1:numel(id)
        x(j)            = proj(st(j), id(j));
        Cost(j)         = maX(st(j));
        nsp(id(j))      = nsp(id(j)) + 1;
       
        % subtract off WtW 
        cproj(j,:) =  proj(st(j) ,:);
        proj(st(j) + trange,:) = proj(st(j) + trange,:)  - x(j) * WtW(:,:,id(j));
    end
    
    xs(i0 + [1:numel(st)])          = x;
    sts(i0 + [1:numel(st)])         = st;
    Costs(i0 + [1:numel(st)])       = Cost;
    ids(i0 + [1:numel(st)])         = id;
    cprojall(i0 + [1:numel(st)], :) = cproj;
    i0 = i0 + numel(st);
end


ids     = ids(1:i0);
xs      = xs(1:i0);
Costs   = Costs(1:i0);
sts     = sts(1:i0);
cprojall = cprojall(1:i0, :);
cprojall = cprojall';

ids = ids - 1;



% keyboard
end %func


%--------------------------------------------------------------------------
function wtw0 =  getWtW2_(Params, W1, W2, utu0)

nt0            = size(W1,1);
W1(2*nt0-1, 1,1) = 0;
W2(2*nt0-1, 1,1) = 0;

fW1 = fft(W1, [], 1);
fW2 = conj(fft(W2, [], 1));

fW2 = permute(fW2, [1 3 2]);

fW = bsxfun(@times, fW1, fW2); 

%
wtw0 = real(ifft(fW, [], 1));

wtw0 = fftshift(wtw0, 1);

utu0 = permute(utu0, [3 1 2]);
wtw0 = bsxfun(@times, wtw0, utu0);
end %func


%--------------------------------------------------------------------------
function Smooth = my_conv_(S1, sig, varargin)

NN = size(S1,1);
NT = size(S1,2);

dt = -4*sig:1:4*sig;
gaus = exp( - dt.^2/(2*sig^2));
gaus = gaus'/sum(gaus);

% Norms = conv(ones(NT,1), gaus, 'same');
%Smooth = zeros(NN, NT);
%for n = 1:NN
%    Smooth(n,:) = (conv(S1(n,:)', gaus, 'same')./Norms)';
%end

Smooth = filter(gaus, 1, [S1' ones(NT,1); zeros(ceil(4*sig), NN+1)]);
Smooth = Smooth(1+ceil(4*sig):end, :);
Smooth = Smooth(:,1:NN) ./ (Smooth(:, NN+1) * ones(1,NN));

Smooth = Smooth';
end %func


%--------------------------------------------------------------------------
function x = sq_(x)
x = squeeze(x);
end %func


%--------------------------------------------------------------------------
function S1 = my_sum_(S1, sig, varargin)
% takes an extra argument which specifies which dimension to filter on
% extra argument can be a vector with all dimensions that need to be
% smoothed, in which case sig can also be a vector of different smoothing
% constants

idims = 2;
if ~isempty(varargin)
    idims = varargin{1};
end
if numel(idims)>1 && numel(sig)>1
    sigall = sig;
else
    sigall = repmat(sig, numel(idims), 1);
end

for i = 1:length(idims)
    sig = sigall(i);
    
    idim = idims(i);
    Nd = ndims(S1);
    
    S1 = permute(S1, [idim 1:idim-1 idim+1:Nd]);
    
    dsnew = size(S1);

    S1 = reshape(S1, size(S1,1), []);
    dsnew2 = size(S1);
    
    S1 = cat(1, 0*ones([sig, dsnew2(2)]),S1, 0*ones([sig, dsnew2(2)]));
    Smax = S1(1:dsnew2(1), :);
    for j = 1:2*sig
        Smax = Smax + S1(j + (1:dsnew2(1)), :);
    end
    
    S1 = reshape(Smax, dsnew);
       
    S1 = permute(S1, [2:idim 1 idim+1:Nd]);
end
end %func


%--------------------------------------------------------------------------
function [score, iY, mu1, mu2, u1, u2] = split_clust_(uu, nhist)

nhist = nhist(:);

nspikes = sum(uu, 1);

uc = zeros(size(uu));
for i = 1:size(uu,2)
    uc(:,i) = my_conv_(uu(:,i)',  max(.5, min(4, 2000/nspikes(i))))'; %.5
%       uc(:,i) = my_conv2(uu(:,i),  max(.25, min(4, 2000/nspikes(i))), 1);
end
%
uc = uc ./repmat(sum(uc,1),size(uc,1), 1);
ucum = cumsum(uc, 1);
%
dd = diff(uc, 1);

iY = zeros(1000,1);
mu1 = zeros(1000,1);
mu2 = zeros(1000,1);
var1 = zeros(1000,1);
var2 = zeros(1000,1);
u1 = zeros(size(uu,1), 1000);
u2 = zeros(size(uu,1), 1000);

maxM = max(uc, [], 1);

inew = 0;

Nfilt = size(uu,2);
mu0 = sum(repmat(nhist(1:100, 1), 1, Nfilt) .* uc, 1);
var0 = sum((repmat(nhist(1:100), 1, Nfilt) - repmat(mu0, 100, 1)).^2 .* uc, 1);

for i = 1:Nfilt
    ix = find(dd(1:end-1, i)<0 & dd(2:end, i)>0);
    
    ix = ix(ucum(ix, i)>.1 & ucum(ix, i)<.8 & uc(ix,i)<.8 * maxM(i)); %.9 not .95
    if nspikes(i) > 500 && numel(ix)>0
        ix = ix(1);
        
        inew = inew + 1;
        
        normuc    = sum(uc(1:ix, i));
        mu1(inew) = sum(nhist(1:ix)     .* uc(1:ix, i))    /normuc;
        mu2(inew) = sum(nhist(1+ix:100) .* uc(1+ix:100, i))/(1-normuc);
        
        var1(inew) = sum((nhist(1:ix)-mu1(inew)).^2     .* uc(1:ix, i))    /normuc;
        var2(inew) = sum((nhist(1+ix:100)-mu2(inew)).^2 .* uc(1+ix:100, i))/(1-normuc);
        
        u1(1:ix,inew) = uu(1:ix, i);
        u2(1+ix:100,inew) = uu(1+ix:100, i);
        
        iY(inew) = i;
    end
    
end

mu1 = mu1(1:inew);
mu2 = mu2(1:inew);
var1 = var1(1:inew);
var2 = var2(1:inew);
u1 = u1(:,1:inew);
u2 = u2(:,1:inew);

n1 = sum(u1,1)';
n2 = sum(u2,1)';
iY = iY(1:inew);

score = 1 - (n1.*var1 + n2.*var2)./((n1+n2).*var0(iY)');
% score = ((n1+n2).*var0(iY)' - (n1.*var1 + n2.*var2))./var0(iY)';
[~, isort] = sort(score, 'descend');

iY = iY(isort);
mu1 = mu1(isort);
mu2 = mu2(isort);
u1 = u1(:,isort);
u2 = u2(:,isort);
score = score(isort);
end %func


%--------------------------------------------------------------------------
function WtW = pairwise_dists_(WU, WUinit)

WU = reshape(WU, [], size(WU,3));
WUinit = reshape(WUinit, [], size(WUinit,3));
WtW = WU' * WUinit;


mu = sum(WU.^2,1);
mu = mu(:);

muinit = sum(WUinit.^2,1);
muinit = muinit(:);

mu     = repmat(mu, 1, size(WUinit,2));
muinit = repmat(muinit', size(WU,2), 1);

WtW = 1 - 2*WtW ./ (muinit + mu);
end %func


%--------------------------------------------------------------------------
function [d2d, iY, drez] = distance_betwxt_(dWU)
[nt0, Nchan, Nfilt] = size(dWU);

dWU = reshape(dWU, nt0*Nchan, Nfilt);
d2d = dWU' * dWU;

mu = sum(dWU.^2,1).^.5;
mu = mu';

muall2 = repmat(mu.^2, 1, Nfilt);
d2d = 1 - 2 * d2d./(1e-30 + muall2+ muall2');

d2d  = 1- triu(1 - d2d, 1);

[dMin, iY] = min(d2d, [], 1);

drez = dMin;
end %func


%--------------------------------------------------------------------------
function v = normc_(v)
v = v./repmat(sum(v.^2, 1), size(v,1),1).^.5;
end %func


%--------------------------------------------------------------------------
function U = zeroOutKcoords_(U, kcoords, criterionNoiseChannels)

[M, imax] = max(abs(U(:,:,1)), [], 1);

% determine over how many channel groups each template exists
aU = sum(U.^2,3).^.5;
ngroups = max(kcoords(:));

aUgroups = zeros(ngroups, size(U,2));
for j = 1:ngroups
    aUgroups(j, :) = mean(aU(kcoords==j,:), 1);
end

% the "effective" number of channel groups is defined below.
% for cases when X channel groups have equal non-zero weights, this number
% equals X
nEffective = sum(aUgroups,1).^2./sum(aUgroups.^2, 1);

[nEffSort, isort] = sort(nEffective, 'descend');

if criterionNoiseChannels<1
    % if this criterion is less than 1, it will be treated as a fraction 
    % of the total number of clusters
    nNoise = ceil(criterionNoiseChannels * size(U,2));
    ThLocal = nEffSort(nNoise);
else
    % if this criterion is larger than 1, it will be treated as the
    % effective number of channel groups at which to set the threshold
    ThLocal = criterionNoiseChannels;
end


for i = 1:size(U,2)
    if ThLocal > nEffective(i)
        U(kcoords~=kcoords(imax(i)),i,:) = 0;
        U(:,i,:) = normc_(squeeze(U(:,i,:)));
    end
end
end %func


%--------------------------------------------------------------------------
function [nS, iNonMatch] = merge_spikes0_(uBase, nS, uS, crit)

if ~isempty(uBase)
    cdot = uBase * uS';
   
    baseNorms = sum(uBase.^2, 2)';
    newNorms  = sum(uS.^2, 2)';
    
    cNorms = 1e-10 + repmat(baseNorms', 1, numel(newNorms)) + repmat(newNorms, numel(baseNorms), 1);
    
    cdot = 1 - 2*cdot./cNorms;
    
    [cdotmin, imin] = min(cdot, [], 1);
    
    iMatch = cdotmin<crit;
    
    nSnew = hist(imin(iMatch), 1:1:size(uBase,1));
    nS = nS + nSnew';
    
    
    iNonMatch = find(cdotmin>crit);
else
   iNonMatch = 1:size(uS,2); 
   nS = [];
end
end %func


%--------------------------------------------------------------------------
function struct_save_(S, vcFile)
save(vcFile, '-struct', 'S', '-v7.3');
fprintf('Wrote to %s\n', vcFile);
end %func


%--------------------------------------------------------------------------
function [uNew, nSnew]= reduce_clusters0_(uS, crit)

 cdot = uS * uS';
 
% compute norms of each spike
newNorms  = sum(uS.^2, 2)';

% compute sum of pairs of norms
cNorms = 1e-10 + repmat(newNorms', 1, numel(newNorms)) +...
    repmat(newNorms, numel(newNorms), 1);

% compute normalized distance between spikes
cdot = 1 - 2*cdot./cNorms;
cdot = cdot + diag(Inf * diag(cdot));

[cmin, newind] = min(single(cdot>crit),[],1);
% if someone else votes you in, your votee doesn't count
% newind(ismember(1:nN, newind)) = [];
newind = unique(newind(cmin<.5));
if ~isempty(newind)
    newind = cat(2, newind, find(cmin>.5));
else
    newind = find(cmin>.5);
end


uNew = uS(newind, :);

nNew = size(uNew,1);

nSnew = merge_spikes0_(uNew, zeros(nNew, 1), uS, crit);
end %func


%--------------------------------------------------------------------------
function [W, U, mu] = get_svds_(dWU, Nrank)

[Wall, Sv, Uall] = svd(gather_try_(dWU), 0);
[~, imax] = max(abs(Wall(:,1)));
Uall(:,1) = -Uall(:,1) * sign(Wall(imax,1));
Wall(:,1) = -Wall(:,1) * sign(Wall(imax,1));

%     [~, imin] = min(diff(Wall(:,1), 1));
%     [~, imin] = min(Wall(:,1));
%     dmax(k) = - (imin- 20);

%     if dmax(k)>0
%         dWU((dmax(k) + 1):nt0, :,k) = dWU(1:nt0-dmax(k),:, k);
%         Wall((dmax(k) + 1):nt0, :)  = Wall(1:nt0-dmax(k),:);
%     else
%         dWU(1:nt0+dmax(k),:, k) = dWU((1-dmax(k)):nt0,:, k);
%         Wall(1:nt0+dmax(k),:) = Wall((1-dmax(k)):nt0,:);
%     end

Wall = Wall * Sv;

Sv = diag(Sv);
mu = sum(Sv(1:Nrank).^2).^.5;
Wall = Wall/mu;

W = Wall(:,1:Nrank);
U = Uall(:,1:Nrank);
end %func


%--------------------------------------------------------------------------
function rez = merge_posthoc2_(rez)
%fracse = 0.1;
mu = rez.mu;
fracse = rez.ops.fracse;

ops = rez.ops;
LAM = ops.lam(3) * (20./mu).^2;
Nfilt = rez.ops.Nfilt;

Cmerge = Inf *ones(Nfilt);
tfi = rez.iNeigh;
tf = rez.cProj;
clusterIDs = rez.st3(:,2);
%
nSpikes = size(rez.st3,1);

fmax = zeros(nSpikes,1, 'single');
pairs = {};
for testID = 1:Nfilt
    spikesTest = clusterIDs==testID;
%     tfnew = bsxfun(@plus, tf(spikesTest, :), LAM(tfi(:, testID))'.*mu(tfi(:, testID))');
%     tf(spikesTest, :) = bsxfun(@rdivide, tfnew, sqrt(1+LAM(tfi(:, testID)))');
    
    pp = tfi(:, testID);
    pp(pp==testID) = [];
    pairs{testID} = pp;
    [~, isame] = min( abs(tfi(:, testID)-testID));
    fmax(spikesTest, 1) = tf(spikesTest, isame);
end


%
inewclust = 0;
clear iMegaC
picked = zeros(Nfilt, 1);
% tic
while 1    
    [maxseed, iseed] = max(rez.nbins(1:Nfilt) .* (1-picked), [], 1);
%     [maxseed, iseed] = max(mu(1:Nfilt) .* (1-picked), [], 1);
    if maxseed<500
        break;
    end
    picked(iseed) = 1;
    % 21, 69,
    %
%     iseed = 410;
    
    run_list = [iseed];
    pair_list = pairs{iseed};
    strun = find(clusterIDs==iseed);
    
    
    while ~isempty(pair_list)
        %
%         picked_pairs = rez.nbins(pair_list);
        
        [mmax, ipair] = max(rez.nbins(pair_list));
        
        
        if mmax<100
            break;
        end
        
        ipair = pair_list(ipair);
        
        %
        imm = ismember(tfi(:, ipair), run_list);
        if sum(imm)
            %
            new_spikes = find(clusterIDs==ipair);
            f1new = max(tf(new_spikes, imm), [], 2);
            
            f2new = fmax(new_spikes);
            
            f1old = fmax(strun);
            f2old = NaN * ones(numel(f1old), 1, 'single');
            i0 = 0;
            for j = 1:length(run_list)
                ifeat = find(tfi(:, run_list(j))==ipair);
                if ~isempty(ifeat)
                    f2old(i0 + (1:rez.nbins(run_list(j))),1) = ...
                        tf(clusterIDs==run_list(j), ifeat);
                    i0 = i0 + rez.nbins(run_list(j));
                end
            end
            
            f1old(isnan(f2old))=[];
            f2old(isnan(f2old))=[];
            mo = merging_score_(f1old - f2old, f1new-f2new, ops.fracse);
            
            
            if mo<3
                strun = cat(1, strun, new_spikes);
                run_list(end+1) = ipair;
                picked(ipair)   = 1;
                if mmax>300
                    pair_list = unique(cat(1, pair_list, pairs{ipair}));
                    pair_list(ismember(pair_list, run_list)) = [];
                end                
            end
        end
        pair_list(pair_list==ipair) = [];
    end
    
    inewclust = inewclust + 1;
    
    iMegaC{inewclust} = run_list;
%     [sum(picked) run_list]
end

% toc
%

iMega = zeros(Nfilt, 1);
for i = 1:length(iMegaC)
   iMega(iMegaC{i}) = iMegaC{i}(1); 
end
rez.iMega = iMega;
rez.iMegaC = iMegaC;


rez.st3(:,5) = iMega(rez.st3(:,2));
end %func


%--------------------------------------------------------------------------
function steps = merging_score_(fold, fnew, fracse)


troughToPeakRatio = 3;

l1 = min(fnew);
l2 = max(fold);

se = (std(fold) + std(fnew))/2;
se25 = fracse * se;
b2 = [0:se25:-l1];
b1 = [0:se25:l2];

hs1 = my_conv_(histc(fold, b1), 1);
hs2 = my_conv_(histc(-fnew, b2), 1);

mmax = min(max(hs1), max(hs2));

m1 = ceil(mean(fold)/se25);
m2 = -ceil(mean(fnew)/se25);

steps = sum(hs1(1:m1)<mmax/troughToPeakRatio) + ...
    sum(hs2(1:m2)<mmax/troughToPeakRatio);
end %func


%--------------------------------------------------------------------------
function [viClu, iMega, iMegaC] = merge_posthoc2_jrc2_(S_clu)
% rez: iNeigh, cProj

fracse = 0.1;
[tfi, tf] = calc_iNeigh_cProj_(S_clu);
nClu = S_clu.nClu;
nbins = histc(S_clu.viClu, .5:1:S_clu.nClu+1); %count histograms
clusterIDs = S_clu.viClu;
nSpikes = numel(clusterIDs);
fmax = zeros(nSpikes,1, 'single');
pairs = cell(1, nClu);
for iClu = 1:nClu
%     spikesTest = clusterIDs==testID;
    viSpk_clu1 = S_clu.cviSpk{iClu};
    pp = tfi(:, iClu);
    pp(pp==iClu) = [];
    pairs{iClu} = pp;
    [~, isame] = min(abs(tfi(:, iClu) - iClu));
    fmax(viSpk_clu1, 1) = tf(viSpk_clu1, isame);
end
inewclust = 0;

% clear iMegaC
picked = zeros(nClu, 1);
iMegaC = {};
while 1    
    [maxseed, iseed] = max(nbins(1:nClu) .* (1-picked), [], 1);
    if maxseed<500
        break;
    end
    picked(iseed) = 1;
    run_list = [iseed];
    pair_list = pairs{iseed};
%     strun = find(clusterIDs==iseed);
    strun = S_clu.cviSpk{iseed};
    
    while ~isempty(pair_list)
        [mmax, ipair] = max(nbins(pair_list));
        if mmax<100, break; end
        
        ipair = pair_list(ipair);
        imm = ismember(tfi(:, ipair), run_list);
        if sum(imm)
%             new_spikes = find(clusterIDs==ipair);
            new_spikes = S_clu.cviSpk{ipair};
            f1new = max(tf(new_spikes, imm), [], 2);            
            f2new = fmax(new_spikes);            
            f1old = fmax(strun);
            f2old = NaN * ones(numel(f1old), 1, 'single');
            i0 = 0;
            for j = 1:length(run_list)
                ifeat = find(tfi(:, run_list(j))==ipair);
                if ~isempty(ifeat)
%                     f2old(i0 + (1:rez.nbins(run_list(j))),1) = tf(clusterIDs==run_list(j), ifeat);
                    f2old(i0 + (1:nbins(run_list(j))),1) = tf(S_clu.cviSpk{run_list(j)}, ifeat);
                    i0 = i0 + nbins(run_list(j));
                end
            end
            
            f1old(isnan(f2old))=[];
            f2old(isnan(f2old))=[];
            mo = merging_score_(f1old - f2old, f1new-f2new, fracse);
            if mo<3
                strun = cat(1, strun, new_spikes);
                run_list(end+1) = ipair;
                picked(ipair) = 1;
                if mmax>300
                    pair_list = unique(cat(1, pair_list, pairs{ipair}));
                    pair_list(ismember(pair_list, run_list)) = [];
                end                
            end
        end
        pair_list(pair_list==ipair) = [];
    end
    
    inewclust = inewclust + 1;
    
    iMegaC{inewclust} = run_list;
%     [sum(picked) run_list]
end

% toc
%

iMega = zeros(nClu, 1);
for i = 1:length(iMegaC)
   iMega(iMegaC{i}) = iMegaC{i}(1); 
end
% rez.iMega = iMega;
% rez.iMegaC = iMegaC;

% rez.st3(:,5) = iMega(rez.st3(:,2)); %post-merge index
viClu = iMega(rez.st3(:,2)); %post-merge index
end %func


%--------------------------------------------------------------------------
function [iNeigh, cProj] = calc_iNeigh_cProj_(S_clu)
% cProj: projections of each detected spike onto the principal components of the channels corresponding to the spike's assigned template. The channel order for each template is available in iNeigh.
% iNeigh: for each template, the channels with largest amplitudes are indexed in order (default 12). This indexing is used to sort coefficients in cProj. Notice this is a fundamentally sparse scheme: only the top channels for each template are stored.
    
% [tfi, tf] = calc_iNeigh_cProj_(S_clu);
% mWtW
nt0 = 61; % number of time points
nNeigh = 16; %in ops.nNeigh
nClu = S_clu.nClu;
Nrank = 3;

error('calc_iNeigh_cProj_: not implemented');

% nspikes = numel(S_clu.viClu);
% nspikes = zeros(nClu, Nbatch);

% cProj = zeros(5e6, nNeigh, 'single');
% nClu = nClu;
% nt0     = rez.ops.nt0;
% Nrank   = ops.Nrank;
WtW     = zeros(nClu,nClu,2*nt0-1, 'single');
for i = 1:Nrank
    for j = 1:Nrank
        utu0 = U0(:,:,i)' * U0(:,:,j);
        if ops.GPU
            wtw0 =  gather_try_(mexWtW2(Params, W(:,:,i), W(:,:,j), utu0));
        else
            wtw0 =  getWtW2_(Params, W(:,:,i), W(:,:,j), utu0);
            wtw0 = permute(wtw0, [2 3 1]);
        end
        WtW = WtW + wtw0;
        clear wtw0 utu0
    end
end
mWtW = max(WtW, [], 3);
WtW = permute(WtW, [3 1 2]);


% sort pairwise templates
% nsp = sum(nspikes,2); %total number of spikes per batch
nsp = numel(S_clu.viClu);
vld = single(nsp>100);
cr    = mWtW .* (vld * vld');
cr(isnan(cr)) = 0;
[~, iNgsort] = sort(cr, 1, 'descend');

% save full similarity score
% rez.simScore = cr;
iNeigh = iNgsort(1:nNeigh, :);

maskTT = zeros(nClu, 'single');
for i = 1:nClu
    maskTT(iNeigh(:,i),i) = 1;
end
end %func
