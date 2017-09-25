
function jrclust(varargin)
% James Jun 2016 01 18 17:05
% addpath('./util');
% [Usage]
% jrclust [command] [parameter_prm.m] -switches
% handles .spk, _evt.mat and _clu.mat files.
% write log to .log file
% as paramters for now. integrate to GUI
%
% [dependencies not yet included]
% jrclust_sort    %GPU
% jrclust_reclust
% plotCluWav
% projectFeatures

% optional switches (at the end)
if nargin==0, dispHelp_(); return; end
if any(ismember(varargin, {'-h', '--help', 'help'})), dispHelp_(); return; end
if any(ismember(varargin, {'-v', '--version'})), dispVer(); return; end
if any(ismember(varargin, {'-p', '--profiler'})), profile off; profile on; end %put it at the end
% if any(ismember(varargin, {'-s', '--save'})), fSaveSpk=1; else fSaveSpk=0; end %cache the result
if any(ismember(varargin, 'commit')), commit_version_(); return; end;
if any(ismember(varargin, 'update')), update_version_(); return; end;
        
% cannot have a switch with -a or -m since I have cluster-auto and
% cluster-manual
% later move to workspace after function exits
Sclu=[]; Sevt=[]; mrWav=[]; mrLfp = []; mrAux = [];
S0 = get(0, 'UserData');
% Subcommand 1/2
vcCommand = varargin{1};
% [csFilter, vcPrompt] = command_filter_(vcCommand);
if nargin >= 2
    vcFile = varargin{2};
elseif nargin == 1
    try
        P = S0.P;
    catch
        P = [];
    end
    if ~isempty(P)
        vcFile = P.vcFile_prm;
        disp(vcFile);
    else
        vcFile = '';
    end
end
% %     if ~isempty(csFilter), vcFile = filename_format(vcFile); end
% else    
%     if ~isempty(vcPrompt)
%         vcFile = uiFileDialog_(csFilter, vcPrompt);
%         if isempty(vcFile), return; end
%     else
%         vcFile = '';
%     end
% end

if numel(varargin) >= 3, vcFile_prb = varargin{3}; else vcFile_prb = ''; end
if numel(varargin) >= 4, vcFile_template = varargin{4}; else vcFile_template = ''; end

% if numel(varargin) >= 3, vcFile_template = varargin{3}; else vcFile_template = ''; end
   

% [vcFile_spk, vcFile_evt, vcFile_clu, vcFile_prm, vcFile_log] = ...
%     subsFileExt(vcFile, '.spk', '_evt.mat', '_clu.mat', '_prm.m', '.log');
% delete('temp_jrclust_*.m');
switch lower(vcCommand)
    case {'quality', 'quality-metric', 'metric'}
        quality_metric_(vcFile); return;
    case 'batch'
        batch_(vcFile, vcFile_prb); return;
    case {'batch-verify', 'batch-validate'}
        batch_verify_(vcFile, vcFile_prb); return;
    case 'verify-detection'
        verify_detection_(vcFile); return;
    case 'featuremap'
        feature_map_(vcFile); return;
    case {'tsf', 'import-tsf'}
        import_tsf_(vcFile); return;
    case {'openephys', 'import-openephys'}
        import_openephys_(vcFile); return;
    case {'abf', 'import-abf'}
        import_abf_(vcFile, vcFile_prb); return;        
    case {'ncs', 'import-ncs', 'import-nlx'}
        import_ncs_(vcFile); return;
    case {'import', 'makegt', 'make-gt'}
        import_clusters_(vcFile); return;
    case 'export'
        export_workspace_(); return;
    case {'doc', 'documentation', 'docu'}
        show_manual_(); return;
    case 'docx'
        open('JRCLUST manual.docx'); return;
    case 'trackset'
        track_set_(vcFile); return;
    case 'install'
        install_jrclust_(); return;
    case 'clear'
        clear global mrWav;
        set(0, 'UserData', []);        
        return;
    case 'edit'
        edit(vcFile); return;
%     case {'trackcentroid', 'track-centroid'}
%         track_centroid_(vcFile); return;
    case 'cabletest'
        cable_test_(vcFile); return;
    case {'traces', 'trace'}
        plotTraces_jrclust_(vcFile); return;
    case 'traces-lfp'
        plotLfp_jrclust_(vcFile); return;  
    case {'lfpcorr', 'lfpcoh'}
        plotLfpCorr_(vcFile); return;
    case 'traces-ap'
        importTraces_ap_(vcFile); return;        
    case 'traces-ref'
        plotRef_jrclust_(vcFile); return;    
    case 'download'
        download_sample_(vcFile); return;
    case 'compile'
        compile_cuda_(); return;
    case 'traces-aux'
        plotAux_jrclust_(vcFiles)
    case 'describe'
        describeCluFile(vcFile); return;
    case 'export-imec-sync'
        export_imec_sync_(vcFile); return;
    case {'cluster-manual', 'manual'}
        clusterManual(vcFile); return;
%     case 'auto-manual'
%         jrclust('auto', vcFile);
%         clusterManual(vcFile); return;
    case 'project'
        projectFeatures(vcFile); return;
    case 'probe'
        plotProbe_(vcFile); 	return;    
    case 'maketrial'
        % ask for channel number and output file
        % confirmation display
        make_trial_(vcFile);
        return;
    case {'makeprm', 'createprm', 'makeprm-f'}
%         P = file2struct('default.prm');  %P = defaultParam();
%         P = appendStruct(P, file2struct('./template.prm')); %overwrite template
        if strcmpi(vcCommand, 'makeprm-f'), fAsk = 0; else fAsk = 1; end
        [P, vcPrompt] = create_prm_file_(vcFile, vcFile_prb, vcFile_template, fAsk);   
        S0.P=P; set(0, 'UserData', S0);
%         disp(vcPrompt);
%         msgbox(vcPrompt, 'modal');
        return;
    case 'thresh'
        Sevt = load_evt_(vcFile);
        figure('color','w', 'Name', vcFile, 'ToolBar', 'none', 'NumberTitle','off'); 
        subplot(211); stem(Sevt.vrThresh_uV); grid on; axis tight;
        xlabel('Site #'); ylabel('Threshold (uV)');
        subplot(212); cdfplot(single(Sevt.vrThresh_uV)); xlabel('Threshold (uV)');
        vrThresh_uV = Sevt.vrThresh_uV;
        assignWorkspace(vrThresh_uV, Sevt);
        return;
    case 'benchtest'
        benchtest_(vcFile); 
        return;
    case {'sitestat', 'sitestats'}
        sitestat_(vcFile);
        return;
    case {'cluster-verify', 'verify', 'validate'}
        clusterVerify_(vcFile);
        return;
    case 'groundtruth'
        groundtruth_(vcFile);
        return;
    case 'import-celldb' %catalin celldb.csv
        import_celldb_(vcFile);        
        return;
    case {'manual-verify'}
        clusterVerify_(vcFile, 1);
        return;        
    case 'siteoffset'
        site_offset_(vcFile); return;
        
    case 'loadraw'
        load_raw_(vcFile); return;
        
    case 'trackdepth'
        track_depth_(vcFile); return;   
        
    case 'event-rate'
        plot_event_rate_(vcFile); return;

    otherwise
        vcFile = subsFileExt(vcFile, '.prm');
        if ~exist(vcFile, 'file'), error('.prm file does not exist'); end
        % load paramter struct
        P = loadParam_(vcFile); %load prm file and create P struct
        if any(ismember(varargin, {'-d', '--debug'})), P.fVerbose=1; end
        if any(ismember(varargin, {'-g1', '--gpu1'})), P.iGpu=1; end
        if any(ismember(varargin, {'-g2', '--gpu2'})), P.iGpu=2; end
        if any(ismember(varargin, {'-g3', '--gpu3'})), P.iGpu=3; end
%         if any(ismember(varargin, {'-s', '--save'})), P.fSaveSpk=1; end %cache the result
        disp('resetting gpu...');
        try gpuDevice(P.iGpu); catch; end
%         gpuDevice();
        
        vcFile_log = subsFileExt(vcFile, '.log');
        
        fid = fopen(vcFile_log, 'w'); 
        warning off;

end

try
% Subcommand 2/2
switch lower(vcCommand)
    case 'activity'
        plot_activity_(P);
        
    case {'raster', 'psth'}
        plot_raster(P);
        
    case 'syncvid'
        led_sync_(P);
        
    case 'showtraj'
        showTraj_(P);
        
    case 'exportcsv'
        % output time(sec), unitid, max chan
        export_csv_(P);
        
    case 'exportcsv0'
        % output time(sec), unitid, max chan
        export_csv_(P, 1);
        
    case 'realtime'
        ui_realtime(P);
        
    case 'filter'
        [mrWav, mrLfp, mrAux, viSiteZero] = export_spk_(P.vcFile, P);
        P.viSiteZero = viSiteZero;

    case 'test' %detect only
        P.fRun=0;
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 0;
        P.fOverwriteFeature = 0;
        P.fOverwriteClu = 0;
        P.fRunCluster = 0;     
        P.fRunMerge = 0;
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);
        
    case 'detectsort' %spike detect and sort       
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 1;
        P.fOverwriteFeature = 1;
        P.fOverwriteClu = 1;
        P.fRunCluster = 1;      
        P.fRunMerge = 1;
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);
        
    case 'detect' %detect only       
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 1;
        P.fOverwriteFeature = 1;
        P.fOverwriteClu = 1;
        P.fRunCluster = 0;        
        P.fRunMerge = 0;
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);
        
    case {'feature'} %feature extraction only       
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 0;
        P.fOverwriteFeature = 1;
        P.fOverwriteClu = 1;
        P.fRunCluster = 0; 
        P.fRunMerge = 0;
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);

%     case 'cluster' %cluster only
%         P.fOverwriteSpk = 0;
%         P.fOverwriteEvt = 0;
%         P.fOverwriteFeature = 0;
%         P.fOverwriteClu = 0;
%         P.fRunCluster = 1;
%         P.fRunMerge = 1;
%         [Sclu, Sevt, mrWav] = process(fid, P);
        
    case {'cluster-auto', 'cluster', 'auto', 'auto-manual', 'auto-verify'} %cluster only
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 0;
        P.fOverwriteFeature = 0;
        P.fOverwriteClu = 1;
        P.fRunCluster = 1;
        P.fRunMerge = 1;
        P.fRunManual = strcmpi(vcCommand, 'auto-manual');
        P.fVerify = strcmpi(vcCommand, 'auto-verify');
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);        
        
    case {'featuresort', 'feature-verify'} %cluster only
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 0;
        P.fOverwriteFeature = 1;
        P.fOverwriteClu = 1;
        P.fRunCluster = 1;
        P.fRunMerge = 1;               
        P.fVerify = strcmpi(vcCommand, 'feature-verify');
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);                
        
    case {'spikesort', 'spikesort-verify'} %detect and cluster
        P.fOverwriteSpk = 1;
        P.fOverwriteEvt = 1;
        P.fOverwriteFeature = 1;
        P.fOverwriteClu = 1;
        P.fRunCluster = 1;
        P.fRunMerge = 1;
        P.fVerify = strcmpi(vcCommand, 'spikesort-verify');
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);        

    case 'cluster-merge'
        P.fOverwriteSpk = 0;
        P.fOverwriteEvt = 0;
        P.fOverwriteFeature = 0;
        P.fOverwriteClu = 1;
        P.fRunCluster = 0;
        P.fRunMerge = 1;
        [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P);

    otherwise
        fprintf(2, 'Invalid command\n');
        dispHelp_(); return; 
end
catch
   disp(lasterr);
    return;
end

% assign to workspace
% assignWorkspace(Sevt, Sclu, mrWav, mrLfp, mrAux);
if isfield(P, 'fRunManual')
    fSet_S0 = ~P.fRunManual;
else
    fSet_S0 = 0;    
end
if fSet_S0
    try
        S0.P=P; S0.Sevt=Sevt; S0.Sclu=Sclu; S0.mrLfp=mrLfp; S0.mrAux=mrAux; 
        set(0, 'UserData', S0);
    catch    
    end
end

% save result, should not duplicate memory
% S0 = get(0, 'UserData');
% S0 = appendStruct(S0, makeStruct(P, Sclu, mrWav, Sevt, mrAux, mrLfp));
% set(0, 'UserData', S0);

% finish up
fclose(fid); %close the log file
disp(['Logged to ', vcFile_log]);
% delete('temp_jrclust_*.m'); %in case not deleted
if ismember(varargin, {'-p', '--profiler'}), profile off; profile viewer; end
end %func


function [csFilter, vcPrompt] = command_filter_(vcCommand)
switch vcCommand
    case 'makeprm'
        csFilter = {'*.bin;*.dat'};
        vcPrompt = 'Select a raw recording ';    
    case 'describe'
        csFilter = {'*.prm'; '*_clu.mat'};
        vcPrompt = 'Select a .prm or clu file';    
    case 'probe'
        csFilter = {'*.prb;'; '*.prm'};
        vcPrompt = 'Select a .prb or .prm file';        
    case {'download', 'cabletest'}
        csFilter = [];
        vcPrompt = '';
    otherwise
        csFilter = {'*.prm'};
        vcPrompt = 'Select a .prm file';
end %swqitch
end %func


function make_trial_(vcFile_prm)
% make a _trial.mat file. stores real time
P = loadParam_(vcFile_prm);

% ask which channel and which output file
csAns = inputdlg('Which channel to load', 'Channel', 1, {num2str(P.nChans)});
iChan = str2double(csAns{1});

% get output file
vcFile_trial = subsFileExt(P.vcFile,  sprintf('_ch%d_trial.mat', iChan));
try
    [FileName,PathName,FilterIndex] = uiputfile(vcFile_trial, 'Save file name');
    if ~FilterIndex, return; end %cancelled
    vcFile_trial = [PathName, FileName];
catch
    fprintf('uiputfile error (old Matlab version). Accepting default');
end

hMsg = msgbox_open('Calculating...');

% load file
fid = memmapfile(P.vcFile, 'Offset', 0, 'Format', P.vcDataType, 'Repeat', inf);
vrWav = fid.Data(iChan:P.nChans:end);
clear fid;

% find rising crossing
if isempty(P.thresh_trial)
    maxV = max(vrWav);
    minV = min(vrWav);
    thresh = (maxV + minV)/2;
else
    thresh = P.thresh_trial;
end
viT = find(vrWav(1:end-1) < thresh & vrWav(2:end) >= thresh) + 1;
nRefrac_trial = round(P.tRefrac_trial * P.sRateHz);
viT = remove_refrac(viT, nRefrac_trial);
vrT = viT / P.sRateHz;

% save
save(vcFile_trial, 'vrT');
disp(['Saved to ', vcFile_trial]);
fprintf('%d events detected\n', numel(vrT));

% edit .prm file
P_trial.vcFile_trial = vcFile_trial;
edit_prm_file_(P_trial, vcFile_prm);

% plot
vlOver = vrWav >= thresh;
figure; hold on; 
stem(vrT, vrWav(viT), 'r');
ylim([minV, maxV]);
plot(find(vlOver)/P.sRateHz, vrWav(vlOver), 'b.');
grid on;
plot(get(gca, 'XLim'), double(thresh) * [1 1], 'r-');
title(sprintf('%d events detected', numel(vrT)));
xlabel('Time (s)'); ylabel(sprintf('Chan %d', iChan)); 
set(gcf, 'Name', vcFile_trial);

msgbox_close(hMsg);

end


function vcFile = uiFileDialog_(vcFilter, vcPrompt)
[vcFile, vcDir, FilterIndex] = uigetfile(vcFilter, vcPrompt);
if FilterIndex==0, vcFile = ''; return ;end
vcFile = [vcDir, vcFile];
disp(vcFile);
end


function Sclu = load_clu_(vcFile_prm)
% Sclu = load_clu_(vcFile_prm)
% Sclu = load_clu_(P)
Sclu = [];
P = loadParam_(vcFile_prm);
try
    if ~matchFileExt(vcFile_prm, '.prm'), error('must provide .prm file'); end
    vcFile_clu = subsFileExt(vcFile_prm, '_clu.mat');
    if ~exist(vcFile_clu, 'file'), vcFile_clu = subsFileExt(P.vcFile, '_clu.mat'); end
    if ~exist(vcFile_clu, 'file'), Sclu=[]; return ;end

    % load from cache
    S0 = get(0, 'UserData');
    if ~isempty(S0)
        if isfield(S0, 'Sclu')
            if strcmpi(S0.Sclu.P.vcFile_prm, vcFile_prm) && ~isempty(S0.Sclu)
                Sclu = S0.Sclu; 
                disp('Sclu is loaded from RAM cache.');
                return; 
            end
        end
    end
    if exist(vcFile_clu, 'file')
        Sclu = load(vcFile_clu); %load Sclu    
        if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end
        disp('Sclu loaded from file.');
        set0(Sclu);
    end
catch
    disp(['load_clu_:', lasterr()]);
end
end %func


function Sgt = load_gt_(vcFile_gt)
% Sgt contains viTime and viClu
if ~exist(vcFile_gt, 'file'), Sgt=[]; return; end
S = load(vcFile_gt);
if isfield(S, 'Sgt')
    Sgt = S.Sgt;
elseif isfield(S, 'viClu') && isfield(S, 'viTime')
    Sgt = S;
elseif isfield(S, 'viClu') && isfield(S, 'viSpk')
    Sgt.viTime = S.viSpk;    
    Sgt.viClu = S.viClu;    
else
    % Convert Nick's format to JRCLUST fomat
    if isfield(S, 'gtTimes')
        Sgt.viTime = cell2mat(S.gtTimes');
        Sgt.viClu = cell2mat(arrayfun(@(i)ones(size(S.gtTimes{i}))*i, 1:numel(S.gtTimes), 'UniformOutput', 0)');
    else
        error('no field found.');
    end
    [Sgt.viTime, ix] = sort(Sgt.viTime, 'ascend');
    Sgt.viClu = Sgt.viClu(ix);
end
[viClu_unique, ~, viClu] = unique(Sgt.viClu);
if max(Sgt.viClu) > numel(viClu_unique)
    Sgt.viClu = viClu;
end
end %func


function groundtruth_(vcFile_prm)
fUseRaw = 0;
P = loadParam_(vcFile_prm);
% P.sRateHz = 20000;
% P.nChans = 60; %assume this
% P.spkLim = [-20, 20];
% P.uV_per_bit = 1;
P.maxAmp = 100;
P.LineStyle = 'k';

if ~isempty(vcFile_prm)
    [vcDir, ~,~] = fileparts(vcFile_prm);
    vcDir = [vcDir, filesep()];
else
    vcDir = '';
end
Sgt = load(subsFileExt(P.vcFile, '_gt.mat'));
if isfield(Sgt,'Sgt'), Sgt=Sgt.Sgt; end
% Sgt = struct();
% vcFile_time = uiFileDialog_([vcDir, '*.csv'], 'PROVIDE CSV or mat FILE');
% if isempty(vcFile_time), return; end
% if matchFileExt(vcFile_time, '.csv')
%     mrCsv = importdata(vcFile_time);
%     Sgt.viTime = round(mrCsv(:,1) * P.sRateHz );
%     Sgt.viClu = mrCsv(:,2);
%     Sgt.viSite = mrCsv(:,3);
% elseif matchFileExt(vcFile_time, '.mat')
%     Sgt = load(vcFile_time);    
%     if isfield(Sgt, 'viSpk'), Sgt.viTime = Sgt.viSpk; end
% else
%     error('not supported');
% end

% vcFile_bin = uiFileDialog_([vcDir, '*.bin'], 'PROVIDE BIN FILE');
% if isempty(vcFile_bin), return; end
if fUseRaw
    P.nSamples_load = 25000*128;
    P.vcCommonRef = 'none';
    % P.fTranspose_bin = 0;
    P.freqLim=[];
    P.nDiff_filt = [];
    P.fft_thresh=0;
    mrWav = import_whisper(P.vcFile, P); %raw channel read
else
    vcFile_spk = subsFileExt(P.vcFile_prm, '.spk'); %file cache
    mrWav = loadWavFile(vcFile_spk, numel(P.viSite2Chan), inf, 'int16');
end
% cvrFet_clu = getFet_clu_(Sgt, mrWav, P); % obtain at max chan. plot cv of fet as a function of SNR
P.nDiff_filt = 0; %do not integrate back
P.fMeanSubt = 0;
[tmrCluWav_mu] = Sclu_meanWav(Sgt, mrWav, P);
P.vcCluWavMode='stdev';
tmrCluWav_sd = Sclu_meanWav(Sgt, mrWav, P);
figure(102); clf; hold on; % show clusters
hPlot = plot_tmr_cluster(tmrCluWav_mu, P);
hold on; plot_tmr_cluster(tmrCluWav_mu+tmrCluWav_sd, P);
hold on; plot_tmr_cluster(tmrCluWav_mu-tmrCluWav_sd, P);
title('black: Gt; Red: Match1; Blue: Match2');
try
    viSite_clu = arrayfun(@(iClu)mode(Sgt.viSite(Sgt.viClu==iClu)), 1:max(Sgt.viClu));
    hold on; plot(viSite_clu, 'r*');
catch
    ;
end
Sgt.tmrWav_clu = tmrCluWav_mu;
Sclu_gt = Sclu_quality(Sgt, P);
% vrVpp_clu = squeeze(max(max(tmrCluWav_mu,[],1)-min(tmrCluWav_mu,[],1),[],2));
% try
%     Sevt = load_evt(vcFile_prm);
%     vnSite_gt = sum(bsxfun(@lt, squeeze(min(tmrCluWav_mu,[],1)), -Sevt.vrThresh_uV),1);
%     vrSnr_gt = vrVpp_clu ./ Sevt.vrThresh_uV();
% catch
%     vnSite_gt = [];
% end
assignWorkspace(Sclu_gt);
% tr1=mr2tr3(mrWav, P.spkLim, Sgt.viTime(Sgt.viClu==3));
% figure; plot(squeeze(single(tr1(:,3,:))))
end %func


function clusterVerify_(vcFile_prm, fReclust)
global mrWav

if nargin<2, fReclust=0; end
fPlotGT = 0;

% if ~matchFileExt(vcFile_prm, '.prm')
%     fprintf(2, 'Must provide a .prm file\n');
%     return;
% end
hMsg=msgbox('verifying...', 'modal');
S0 = get(0, 'UserData');
try
    P = loadParam_(vcFile_prm);    
    mrWav = load_spk_(vcFile_prm);
    Sclu = load_clu_(vcFile_prm);
    if fReclust
        Sclu = post_merge(Sclu, mrWav, P);
    end
        
    vcFile_gt = subsFileExt(vcFile_prm, '_gt.mat'); %ground truth file    
    if ~exist(vcFile_gt, 'file')
        vcFile_gt = subsFileExt(P.vcFile, '_gt.mat');
    end
    Sgt = load_gt_(vcFile_gt);    
    if isempty(Sgt), fprintf(2, 'Groundtruth does not exist. Run "jrclust import" to create a groundtruth file.\n'); return; end
    disp('Using manually curated clusters');
    
    % plot GT
    if fPlotGT
        tmrWav_clu = meanCluWav(mrWav, Sgt.viClu, Sgt.viTime, P);
        P.maxAmp = 1;
        figure; plot_tmr_cluster(tmrWav_clu, P); %do this after plotSpk_
    end
    
    if isempty(Sclu), Sclu = Sgt; disp('Self-check, expecting 100% accuracy'); end
%     Sclu.viClu(isnan(Sclu.delta)) = 0; %debug. will minimize false positives
    fprintf('verifying cluster...\n'); 
    [mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = ...
        clusterVerify(Sgt.viClu, Sgt.viTime, Sclu.viClu, Sclu.viTime, 20);  %Sgt.viTime
    
    vcFile_score = strrep(vcFile_prm, '.prm', '_score.mat');
    S_score = makeStruct(mrMiss, mrFp, vnCluGt, miCluMatch, P, Sgt, S_score_clu);
    S_score.cviTime_clu = Sclu.cviSpk_clu(S_score_clu.viCluMatch)';
    S_score.cviTime_gt = arrayfun(@(iClu)Sgt.viTime(Sgt.viClu==iClu), 1:max(Sgt.viClu), 'UniformOutput', 0);
    
   
    % display groundtruth Vpp and SNR
    P.nDiff_filt=0;
    trWav_gt = Sclu_meanWav(Sgt, mrWav, P);
    [S_score.vrVpp_gt, viSite_gt] = max(shiftdim(max(trWav_gt,[],1)-min(trWav_gt,[],1)),[],1);
    [S_score.vrVmin_gt, ~] = min(shiftdim(min(trWav_gt,[],1)));
    Sevt = load_evt(P.vcFile_prm);
    S_score.vrVrms_site = Sevt.vrThresh_uV / Sevt.P.qqFactor;
    vr_Vrms_site_gt = S_score.vrVrms_site(viSite_gt(:));
    S_score.vrSnr_gt = abs(S_score.vrVpp_gt(:)) ./ vr_Vrms_site_gt;
    S_score.vrSnr_min_gt = abs(S_score.vrVmin_gt(:)) ./ vr_Vrms_site_gt;
    
    S_score.trWav_gt = trWav_gt;
    S_score.viSite_gt = viSite_gt;
    S_score.S_metric = quality_metric_(Sevt, Sclu, mrWav);
    S_score.vrSnr_min_clu = S_score.S_metric.vrSnr_clu(S_score_clu.viCluMatch);
    
% %     S_score.vrFetCv_clu = Sclu_FetCv(S_score_clu.cviHit_clu, S_score.viSite_gt, mrWav, P); % quantify feature stability by GT unit
    mrAmin_gt = shiftdim(min(S_score.trWav_gt,[],1));
    S_score.vnSite_gt = sum(bsxfun(@lt, mrAmin_gt, -Sevt.vrThresh_uV),1);
    fprintf('SNR_gt (Vpp/Vrms): %s\n', sprintf('%0.1f ', S_score.vrSnr_gt));
    fprintf('SNR_gt (Vp/Vrms): %s\n', sprintf('%0.1f ', S_score.vrSnr_min_gt));
    fprintf('nSites>thresh (GT): %s\n', sprintf('%d ', S_score.vnSite_gt));
%     
    save(vcFile_score, '-struct', 'S_score');
    fprintf('Saved to %s.\n', vcFile_score); 

    %     figure; plotRD_cluster_(Sclu, miCluMatch(1,:)); %best matching cluster
%     fprintf('took %0.1fs\n', toc(t2));

    try close(hMsg); catch; end
%     return;

%     vrVpp_clu = max(Sclu.tmrWav_clu)-min(Sclu.tmrWav_clu);
%     vrVpp_clu = max(shiftdim(vrVpp_clu), [], 1);    
%     vrVpp_clu = vrVpp_clu(S_score_clu.viCluMatch);
    set0(S_score);       
    assignWorkspace(S_score); %put in workspace
    
    figure; set(gcf,'Name',P.vcFile_prm);
    subplot 121; plot_cdf(S_score.S_score_clu.vrFp); hold on; plot_cdf(S_score.S_score_clu.vrMiss); 
    legend({'false positives', 'miss rates'}); ylabel('CDF'); grid on; xlabel('Cluster count');
    
    subplot 122; hold on;
    plot(S_score.vrSnr_min_gt, S_score.S_score_clu.vrFp, 'b.', S_score.vrSnr_min_gt, S_score.S_score_clu.vrMiss, 'r.');
    legend({'false positives', 'miss rates'}); ylabel('score'); grid on; xlabel('SNR (Vp/Vrms)');
    
%     subplot 133; hold on;
    disp('SNR>8 stats');
    viSnr8_gt = find(S_score.vrSnr_min_gt>8);
    fprintf('\tfalse positive (%%): '); disp_stats(S_score.S_score_clu.vrFp(viSnr8_gt)*100);
    fprintf('\tfalse negative (%%): '); disp_stats(S_score.S_score_clu.vrMiss(viSnr8_gt)*100);
%     legend({'false positives', 'false negative'}); ylabel('score'); grid on; xlabel('SNR (Vp/Vrms)');
%     figure; multiplot([],300,[],trWav_gt); title('GT waveform (scale:300)');
    return;
    
%     
%     if ~isfield(Sgt, 'vrVpp_clu'), return ;end
% %     Sgt.tmrWav_clu = meanCluWav(mrWav, Sgt.viClu, Sgt.viTime, P);
%     S_score.vrVpp_clu = Sgt.vrVpp_clu;
% %     S_score.vrVpp_clu
%     
%     set0(S_score);    
%     vl100uV = S_score.vrVpp_clu>100;
%     fprintf('(25-50-75%%)>100uVpp: <score>:%0.1f %0.1f %0.1f; <miss>:%0.1f %0.1f %0.1f; <FP>:%0.1f %0.1f %0.1f\n', ...
%         quantile(S_score_clu.vrScore(vl100uV),[.25,.5,.75])*100, ...
%         quantile(S_score_clu.vrMiss(vl100uV),[.25,.5,.75])*100, ...
%         quantile(S_score_clu.vrFp(vl100uV),[.25,.5,.75])*100);

    
    try
    figure; 
    subplot 221; plot(S_score.vrVpp_clu, S_score_clu.vrScore, '*'); 
    ylabel('1-FP-missed'); grid on; xlabel('Vpp amplitude of matching cluster'); 
    title(P.vcFile);

    subplot 223; plot(S_score.vrVpp_clu, 1-S_score_clu.vrFp, '*'); 
    ylabel('Purity'); grid on; xlabel('Vpp amplitude of matching cluster'); 
    subplot 224; plot(S_score.vrVpp_clu, 1-S_score_clu.vrMiss, '*'); 
    ylabel('Completeness'); grid on; xlabel('Vpp amplitude of matching cluster'); 
    subplot 222; plot(1-S_score_clu.vrFp, 1-S_score_clu.vrMiss, '*'); 
    xlabel('Purity'); ylabel('Completeness'); grid on;
    set(gcf,'Color','w', 'Name', P.vcFile);  
    catch
        disp(lasterr());
    end
%     figure; S_score_plot(S_score);
    return;

    mrWav = import_whisper(P.vcFile, 'nChans', P.nChans); %raw channel read
    figure(102); clf; hold on; % show clusters
%     mrWav = loadWav_prm_(vcFile_prm);
%     mrWav = get0('mrWav');
    P.uV_per_bit = 100;
    tmrCluWav_mu = meanCluWav(mrWav, Sgt.viClu, Sgt.viTime, P);
%     tmrCluWav_mu = permute(mr2tr3(mrWav, P.spkLim, S.gtTimes{2}),[1,3,2]); %Sgt.viTime(Sgt.viClu==3)
    P.LineStyle = 'k';
    hPlot = plot_tmr_cluster(tmrCluWav_mu, P);
    title('black: Gt; Red: Match1; Blue: Match2');
%     Sclu1 = struct('cl', Sgt.viClu, 'viTime', Sgt.viTime, 'P', Sclu.P);
%     uiCluWav(mrWav, Sclu1, P);
    
    % matching clu. only plot matching clu
    iClu1 = miCluMatch(1,1);
    iClu2 = miCluMatch(2,1);   
    viPlot = ismember(Sclu.viClu, iClu1);
    tmrCluWav_mu1 = meanCluWav(mrWav, Sclu.viClu(viPlot), Sclu.viTime(viPlot), P);
    P.LineStyle = 'r';
    hPlot1 = plot_tmr_cluster(tmrCluWav_mu1, P);
    

    % second closest match
    viPlot = ismember(Sclu.viClu, iClu2);
    tmrCluWav_mu1 = meanCluWav(mrWav, Sclu.viClu(viPlot), Sclu.viTime(viPlot), P);
    P.LineStyle = 'b';
    hPlot1 = plot_tmr_cluster(tmrCluWav_mu1, P);

     % projection
    figure(103); clf;
%     Sclu.viClu = Sclu.cl; 
    P.fMeanSiteRef = 0;
    mhAx = plotProj(Sclu, mrWav, P, iClu1, iClu2);

catch
    disp(lasterror());
end
try close(hMsg); catch; end
end %func


function plotRD_cluster_(Sclu, viCluMatch)
viClu = Sclu.icl(viCluMatch);
vrX = double(log10(Sclu.rho(viClu)));
vrY = double(log10(Sclu.delta(viClu)));
hold on; plot(vrX, vrY, 'ro');
csClu = arrayfun(@(i)sprintf('%d', i), 1:numel(viClu), 'UniformOutput', 0);
text(vrX, vrY, csClu);
end


function export_csv_(P, fZeroIndex)
if nargin<2, fZeroIndex = 0; end
% output time(sec), unitid, max chan
% if ~matchFileExt(vcFile_prm, '.prm')
%     fprintf(2, 'Must provide a .prm file\n');
%     return;
% end

% P = loadParam(vcFile_prm);
vcFile_clu = subsFileExt(P.vcFile_prm, '_clu.mat');
Sclu = load(vcFile_clu); %load Sclu    
if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end

vrTime = double(Sclu.viTime) / Sclu.P.sRateHz;
viClu = double(Sclu.viClu);
viSite = double(Sclu.viSite) - fZeroIndex; %zero base

vcFile_csv = subsFileExt(P.vcFile_prm, '.csv');
dlmwrite(vcFile_csv, [vrTime(:), viClu(:), viSite(:)], 'precision', 9);
% csvwrite(vcFile_csv, [vrTime(:), viClu(:)]);
fprintf('wrote to %s\n', vcFile_csv);
end %func


function dispHelp_()
    csHelp = {...
        'JRCLUST (last updated 2016 July 22)'; %@TODO: script update
        '   Created by James Jun'; 
        '   Applied Physics and Instrumentation Group';
        '   HHMI - Janelia Research Campus';
        ''; 
        'Usage: jrclust command [filename]';
        '';
        'optional arguments: (attach at the end)';
        '  -h, --help            show this help message and exit ';
        '  --version, -v         print the version of jrclust ';
        '  --debug, -d           activate debug logging mode (default: False) ';
        '  --profiler, -p        activate the profiler (default: False) ';
        '  --gpu#, -g#           Use GPU # (e.g. -g1) ';
        '  --save, -s            save .spk file (cached) ';
        '';
        'Main commands:';
        '  jrclust install';
        '                        Install jrclust by compiling codes';
        '  jrclust clear';
        '                        Clear cache';
        '  jrclust update';
        '                        Update code from the dropbox location'; 
        '  jrclust compile';
        '                        Recompile CUDA code (GPU codes, *.cu)'; 
        '  jrclust doc';
        '                        Open JRCLUST documentation'; 
        '  jrclust download sample';
        '                        Download sample data from Neuropix phase 2 probe';
        '  jrclust probe myprobe.prb or myparams.prm';
        '                        Plot probe layout'; 
        '  jrclust makeprm myrecording.bin myprobe.prb template.prm';
        '                        create a new parameter file based on the template file and probe file';
        '  jrclust makeprm myrecording.bin myprobe.prb (postfix)';
        '                        Generates myrecording(postfix).prm file from a template (default.prm + template.prm). Optionally attach postfix for multi-shank file';
        '  jrclust edit myparams.prm';
        '                        Edit myparams.prm file'                
        '  jrclust traces myparams.prm';
        '                        Displays raw trace';                
        '  jrclust spikesort myparams.prm';
        '                        Run the whole suite (spike detection and clustering) ';
        '  jrclust describe myparams.prm';
        '                        Display information about a clu dataset ';
        '  jrclust manual myparams.prm';
        '                        Run the manual clustering GUI ';
        '  jrclust auto-manual myparams.prm';
        '                        Run the auto clustering and do the manual clustering next';        
        '  jrclust exportcsv myparams.prm';
        '                        Export clustered information to a csv file (spike time, cluster #, max site#)';
        '  jrclust export-imec-sync myparams.prm';
        '                        Export Sync channel (uint16) to the workspace (vnSync)';        
        '';
        'Advanced commands:';
        '  jrclust filter myparams.prm'; 
        '                        Filters the raw recording (generates .spk)'; 
        '  jrclust detect myparams.prm';
        '                        Run spike detection and feature detection (generates _evt.mat)';      
        '  jrclust feature myparams.prm';
        '                        Run feature detection (generates _evt.mat)';                                  
        '  jrclust detectsort myparams.prm';
        '                        Re-run spike detection, feature extraction and sorting';
        '  jrclust featuresort myparams.prm';
        '                        Re-run feature extraction and sorting (generates _evt.mat, _clu.mat)';
        '  jrclust featuresort myparams.prm';
        '                        Run feature detection and sort';     
        '  jrclust cluster myparams.prm';
        '                        cluster the data (after spike detection) ';
        '  jrclust thresh myparams.prm';
        '                        plot spike detection threshold';
        '  jrclust sitestat myparams.prm';
        '                        Site by site analysis';
        '  jrclust siteoffset myparams.prm';
        '                        show the offset of the raw traces'
        '  jrclust cabletest myrecording.bin';
        '                        Phase 3 cable test'
        '  jrclust trackset mycollection.set';
        '                        Analyze set of files'; 
        '  jrclust activity myparams.prm';
        '                        Show firing rate as a function of time and depth';         
        '';
        'Experimental commands:';
        '  jrclust verify myparams.prm';
        '                        Compares against ground truth file (_gt.mat)';     
        '  jrclust trackdepth myparams.prm';
        '                        LFP based depth tracking'        
        '  jrclust psth myparams.prm'; 
        '                        plot peristimulus time histograms for each cluster (trial-averaged)';
        '  jrclust raster myparams.prm'; 
        '                        plot rastergrams for each cluster';    
        '  jrclust syncvid myparams.prm';
        '                        Synchronize video using LED blinking'; 
        '  jrclust project myfile_evt.mat';
        '                        show max-amp projection view';        
        '  jrclust cluster-merge myparams.prm';
        '                        merges clusters based on waveforms ';   
        };
    cellfun(@(s)fprintf('%s\n',s), csHelp);
    
    %show_manual_();
end %func


function show_manual_(vcFile_doc)
if nargin<1, vcFile_doc = 'JRCLUST manual.pdf'; end
S_cfg = read_cfg_();
open([S_cfg.path_dropbox, filesep(), vcFile_doc]);
end




function dispVer()
% version history. maintain version history
disp('0.1.0 (2016 Jan 18)');
end %func


% function downloadData(csArg)
% % download a sample data and quit. save to a current directory
% % connect to dropbox link
% disp(csArg)
% disp('to be implemented');
% end %func
% 
% 
% function plot_psth_(vcFile_prm)
% 
% % plot psth
% P = loadParam(vcFile_prm);
% if isfield(P, 'vcFile_psth'), P.vcFile_trial = P.vcFile_psth; end
%  
% vrTime_trial = loadTrial(P.vcFile_trial);
% 
% % import cluster time
% vcFile_clu = subsFileExt(P.vcFile, '_clu.mat');
% Sclu = loadClu(vcFile_clu);
% viClu = Sclu.viClu;
% vrTime = double(Sclu.viTime) / Sclu.P.sRateHz;
% 
% % plot lcuster times
% cvrCluTime = arrayfun(@(iClu)vrTime(viClu==iClu), 1:max(viClu), 'UniformOutput', 0);
% 
% nClu = numel(cvrCluTime);
% tbin = P.tbin_psth;
% nSmooth = round(P.nSmooth_ms_psth);
% tlim = P.tlim_psth;
% nlim = round(tlim/tbin);
% vrTimePlot = (nlim(1):nlim(end))*tbin;
% % viTrial = 1:numel(vrTime_trial);
% viTime_Trial1 = round(vrTime_trial / tbin);
% if P.fAverageTrial_psth, figure('Name', P.vcFile); hold on; end
% 
% hFig = figure;
% hTabGroup = uitabgroup(hFig);
% 
% for iClu=1:nClu
%     htab1 = uitab(hTabGroup, 'Title', sprintf('Clu %d', iClu));
%     hax1 = axes('Parent', htab1);
%     axes(hax1);
%     
%     vlTime1=zeros(0);
%     vlTime1(ceil(cvrCluTime{iClu}/tbin))=1;
%     if ~P.fAverageTrial_psth
%         figure;
%         mrRate1 = vr2mr2(double(vlTime1), viTime_Trial1, nlim);
%         mrRate1 = filtfilt(ones(nSmooth,1),nSmooth,mrRate1)/tbin;
%         plot(vrTimePlot, mrRate1);
%         title(sprintf('Cluster %d', iClu));         
%     else
%         mr1 = vr2mr2(double(vlTime1), viTime_Trial1, nlim);
%         vrRate1 = mean(mr1,2);
%         vrRate1 = filtfilt(ones(nSmooth,1),nSmooth,vrRate1)/tbin;
%         plot(vrTimePlot, vrRate1);
%     end
%     if ~isempty(P.rateLim_psth), set(gca, 'YLim', P.rateLim_psth); end
%     set(gcf, 'Name', P.vcFile);
%     set(gca,'XTick', tlim(1):P.xtick_psth:tlim(end));   grid on;
%     xlabel('Time (s)'); ylabel('Firing Rate (Hz)'); xlim(tlim);
% end
% if P.fAverageTrial_psth
%     legend(arrayfun(@(i)sprintf('Clu%d', i),nClu:-1:1, 'UniformOutput',0))
% end
% 
% end %func


function [flag, P] = check_cache_(vcFile_prm)
% returns 1 if the cache is valid and current
% load spk from a file
flag = 0; P = [];
try
    if ischar(vcFile_prm)
        if ~matchFileExt(vcFile_prm, '.prm'), error('Must provide .prm file'); end
        P = loadParam_(vcFile_prm);
    elseif isstruct(vcFile_prm)
        P = vcFile_prm;
        vcFile_prm = P.vcFile_prm;
    else
        error('check_cache_: invalid input');
    end

    S0 = get(0, 'UserData');
    if isfield(S0, 'P')
        P0 = S0.P;
    elseif isfield(S0, 'Sevt')
        P0 = S0.Sevt.P;
    elseif isfield(S0, 'Sclu')
        P0 = S0.Sclu.P;
    else
        P0 = [];
    end

    flag = strcmpi(P0.vcFile, P.vcFile);
catch
    flag = 0;
end
end %end


function mrWav = load_spk_(vcFile_prm)
global mrWav
if isstruct(vcFile_prm), vcFile_prm=vcFile_prm.vcFile_prm; end
[fCache, P] = check_cache_(vcFile_prm);
if fCache && ~isempty(mrWav), disp('loaded from cache'); return; end

try
    vcFile_spk = subsFileExt(P.vcFile_prm, '.spk'); %file cache
    if ~exist(vcFile_spk, 'file'), vcFile_spk = replacePath(vcFile_spk, vcFile_prm); end
    if ~exist(vcFile_spk, 'file'), vcFile_spk = subsFileExt(P.vcFile, '.spk'); end
    mrWav = loadWavFile(vcFile_spk, numel(P.viSite2Chan), inf, 'int16');
    set0(P);
%     set0(mrWav); %cache this result
    disp('Loaded from Disk and cached in RAM.');
catch
    disp(lasterr());
    mrWav = [];
end
end %func


function [Sclu, mrWav] = clusterManual(vcFile_prm)
global mrWav

if ~exist(vcFile_prm, 'file')
    fprintf(2, 'file not found. cluster the data first\n');
    return;
end
if ~matchFileExt(vcFile_prm, '.prm')
   fprintf(2, 'Provide a .prm file\n');
   return;
end 
    
P = loadParam_(vcFile_prm);
if ~check_cache_(P)
    mrWav = [];
    set(0, 'UserData', []);
end
vcFile_clu = subsFileExt(P.vcFile_prm, '_clu.mat');
edit(vcFile_prm);

% Load clu file
% Sclu = load(vcFile_clu);
Sclu = load_clu_(vcFile_prm);
if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end
P_clu = Sclu.P;

Sevt = load_evt_(vcFile_prm);
mrWav = load_spk_(vcFile_prm);

% add common ref for display
if P.fAddCommonRef, mrWav = addCommonRef_(mrWav, Sclu.P); end

% merge clu if parameter values changed
% if isempty(P.maxWavCor) ~= isempty(P_clu.maxWavCor)
% if P.rho_cut ~= P_clu.rho_cut || P.delta1_cut ~= P_clu.delta1_cut || P.min_count ~= P_clu.min_count
switch lower(questdlg('Load last saved?', 'Confirmation'))
    case 'no'
        Sclu = post_merge(Sclu, mrWav, P);
    case 'cancel'
        return;
end

% visualize cluster waveforms
set0(Sclu, Sevt, P);
jrclust_ui(mrWav, Sclu, P);
end %func


function mrWav = addCommonRef_(mrWav, P)
switch lower(P.vcCommonRef)
    case {'none', 'holtzman'}    
        return; 
end

fprintf('Adding common reference...'); t1 = tic;
fid = fopen(subsFileExt(P.vcFile, '.ref'));
vrWav_ref = fread(fid, Inf, 'int16=>int16');
fclose(fid);
for iCh=1:size(mrWav,2)
    if ismember(iCh, P.viSiteZero), continue; end
    mrWav(:,iCh) = mrWav(:,iCh) + vrWav_ref;  %memory efficient 
end
% mrWav(:,viSite) = bsxfun(@plus, mrWav(:,viSite), vrWav_ref);
fprintf('took %0.1fs\n', toc(t1));
end %fund


function vcFile1 = fullpath2fname(vcFile)
[~,vcFile1, vcExt1] = fileparts(vcFile); 
vcFile1 = [vcFile1, vcExt1];
end


function [Sclu, Sevt, mrWav, mrLfp, mrAux] = process_(fid, P)
global mrWav

Sclu=[]; Sevt=[]; mrLfp = []; mrAux = []; %deal did not work for some reason
S0 = get(0, 'UserData');
if ~check_cache_(P)
    mrWav=[]; 
    set(0, 'UserData', []); 
    S0 = [];
end
% try
%     if ~strcmpi(S0.P.vcFile, P.vcFile)
%         S0=[];         
%     end
% catch
%     S0=[]; 
% end
% if isempty(S0), set(0, 'UserData', []); end

try    
    ffprintf([1,fid], ['Loading ' P.vcFile]);
    if ~P.fRun, return; end %just list files to be processed   
    vcFile_spk = subsFileExt(P.vcFile_prm, '.spk');
    if ~exist(vcFile_spk, 'file')
        vcFile_spk1 = subsFileExt(P.vcFile, '.spk');
        if exist(vcFile_spk1, 'file'), vcFile_spk = vcFile_spk1; end
    end
    [vcFile_evt, vcFile_clu] = subsFileExt(P.vcFile_prm, '_evt.mat', '_clu.mat');

    if ~exist(vcFile_clu, 'file') || P.fOverwriteClu
        if ~exist(vcFile_evt, 'file') || P.fOverwriteEvt
            fRunEvt = 1;
        else
            try 
                Sevt = S0.Sevt;
            catch
                Sevt = load_evt_(P.vcFile_prm);
            end
            fRunEvt = isempty(Sevt);
        end
        if fRunEvt %event detection
            if ~exist(vcFile_spk, 'file') || P.fOverwriteSpk
                clear global mrWav %memory clear
                [mrWav, mrLfp, mrAux, viSiteZero] = export_spk_(P.vcFile, P);
%                 write_bin_mrWav_(subsFileExt(P.vcFile, '.spk')); %timer write
                P.viSiteZero = viSiteZero;
%                 set_S0(mrWav, mrLfp, mrAux, P);
                ffprintf([1,fid], ['Wrote to ', vcFile_spk]);
            else
                mrWav = load_spk_(P);
            end
            Sevt = jrclust_event_(mrWav, P);
%             set_S0(Sevt);
            if P.fSaveEvt
                save_evt_(vcFile_evt, Sevt);
                ffprintf([1,fid], ['Wrote to ', vcFile_evt]); 
            end                
        else
            if P.fOverwriteFeature
                mrWav = load_spk_(P);
                Sevt = getFet_(Sevt, mrWav, P);                    
                if P.fSaveEvt
                    save_evt_(vcFile_evt, Sevt);                        
                    ffprintf([1,fid], ['Wrote to ', vcFile_evt]); 
                end
            end
        end

        %run cluster        
        if isempty(mrWav), mrWav = load_spk_(P); end
        S0.P=P; S0.Sevt=Sevt;
        set(0, 'UserData', S0);
        if P.fRunCluster
            Sclu = jrclust_sort(Sevt, P);   
        elseif P.fRunMerge %prep merge
            Sclu = loadClu(vcFile_clu);                
        end
        if P.fRunMerge %merge                
            Sclu = post_merge(Sclu, mrWav, P);
        end
        if P.fRunMerge || P.fRunCluster %save
            save_struct_(vcFile_clu, Sclu);
            ffprintf([1,fid], ['wrote to ', vcFile_clu]);
            describeCluFile(Sclu);
        end
        if isempty(Sclu)
            Sclu = load(vcFile_clu);
            if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end
        end
        describeCluFile(Sclu);
        S0.Sclu = Sclu;
        set(0, 'UserData', S0);
        if struct_flag(P, 'fRunManual')            
            jrclust_ui(mrWav, Sclu, P);
        end
        if struct_flag(P, 'fVerify')
            clusterVerify_(P.vcFile_prm);
        end        
    end
catch
    disp(lasterror());
%         ffprintf([2,fid], lasterror());
    disp('resetting gpu...');
    try gpuDevice(P.iGpu); catch; end        
%         gpuDevice();          
end %try
% end %for
end %func


function flag = struct_flag(P, vcFlag)
if ~isfield(P, vcFlag), flag = 0; return; end
flag = P.(vcFlag);
end %func


function corr = wave_corr_(mr1, mr2)
% both are mean subtracted
corr = mean(mr1(:) .* mr2(:)) / std(mr1(:)) / std(mr2(:));
end %func


function vrRef = load_ref_(vcFile_ref, P)
% if second argument specified return index containing valid
tbin_ref = .01; %10 msec bin
try
    fid = fopen(subsFileExt(vcFile_ref, '.ref'), 'r');
    vrRef = fread(fid, inf, 'int16=>single');
    fclose(fid);
    if nargin >= 2
        nwin = round(P.sRateHz * tbin_ref);
        vrRef_bin = std(reshape_vr2mr(vrRef, nwin), 1,1);
        vlRef_bin = thresh_mad(vrRef_bin, P.rejectSpk_mean_thresh); %keep index
        vrRef = repmat(vlRef_bin(:)', [nwin, 1]); 
        vrRef = vrRef(:);

%         vrRef = abs(vrRef);
%         vrRef = vrRef < median(subsample_vr(vrRef, 300000)) * P.reject_mean_thresh;        
%         vrRef = abs(vrRef) < std(vrRef) * reject_mean_thresh;
        fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vrRef))*100 );
    end %if
catch
    disp(['load_ref_: ', lasterr]);
    vrRef = [];
end

end %func


function Sevt = getFet_(Sevt, mrWav, P)
if P.fRejectSpk_vpp || strcmpi(P.vcFet, 'vpp') || strcmpi(P.vcFet, 'energy')
    [mrVpp, mrEnergy] = getVpp(mrWav, Sevt.viSpk, P);
end 
if P.fRejectSpk_vpp    
    % remove spikes in Sevt
    viRemove = find(median(mrVpp) >= mean(Sevt.vrThresh_uV * 2));
    fprintf('Removed %d spikes (%0.1f%%)\n', numel(viRemove), numel(viRemove)/size(mrVpp,2)*100);
    Sevt.viSpk(viRemove) = [];
    Sevt.vrSpk(viRemove) = [];
    Sevt.viSite(viRemove) = [];
    mrVpp(:,viRemove) = [];
    try mrEnergy(:,viRemove) = []; catch; end
elseif ~strcmpi(P.vcFet, 'vpp')
    mrVpp = []; %clear memory if not needed
end
switch lower(P.vcFet)
    case 'pca'
        if ~P.fPcaDetect
            [Sevt.trPv, Sevt.trFet, Sevt.mrPd] = ...
                getPrinVec_gpu(mrWav, Sevt.viSpk, Sevt.viSite, P); %4 uses matlab pca
        else
%             Sevt.trFet = []; %cluster algorithm will project directly
            [Sevt.trFet, Sevt.miSites_fet] = pca_project_(mrWav, Sevt, P); % P.mrPv used
%             Sevt.trFet = pca_project_(mrWav, Sevt, P); % P.mrPv used
        end
    case {'xcor', 'xcov'}
        if ~P.fPcaDetect, error('fPcaDetect must be enabled when using xcor feature'); end
        [Sevt.trFet, Sevt.miSites_fet] = xcor_project_(mrWav, Sevt, P); % P.mrPv used
    case 'moment'
        if ~P.fPcaDetect, error('fPcaDetect must be enabled when using moment feature'); end
        [Sevt.trFet, Sevt.miSites_fet] = moment_project_(mrWav, Sevt, P); % P.mrPv used
    case {'autocor', 'acor'}
%         if ~P.fPcaDetect, error('fPcaDetect must be enabled when using moment feature'); end
        [Sevt.trFet, Sevt.miSites_fet] = autocor_project_(mrWav, Sevt, P); % P.mrPv used
    case 'spacetime'
        [Sevt.trFet, Sevt.miSites_fet] = fet_spacetime_(mrWav, Sevt, P); % P.mrPv used
    case 'tetrode'
        Sevt.trFet = fet_tetrode_(mrWav, Sevt, P); % P.mrPv used
    case 'diff'
        if ~P.fPcaDetect, error('fPcaDetect must be enabled when using diff feature'); end
        [Sevt.trFet, Sevt.miSites_fet] = diff_project_(mrWav, Sevt, P); % P.mrPv used
    case 'spca'
        [Sevt.trFet, Sevt.cviFet_site, Sevt.miSites_site] = ...
            pca_project_site_(mrWav, Sevt, P); % P.mrPv used
    case {'amp', 'amp2'}
        Sevt.trFet = shiftdim(getAmpMin(mrWav, Sevt.viSpk, P), -1);    
        if P.vcFet(end)=='2', Sevt.trFet=Sevt.trFet.^2; end
    case {'vpp', 'vppsqrt'}
        if strcmpi(P.vcFet, 'vppsqrt'), mrVpp = sqrt(abs(mrVpp)); end
        Sevt.trFet = shiftdim(mrVpp, -1);
    case 'energy'
        Sevt.trFet = shiftdim(mrEnergy, -1);
    case {'ampnorm', 'ampnorm2'} %normalize by threhold
        mrAmp = getAmpMin(mrWav, Sevt.viSpk, P);
%         vrSd = std(single(mrWav(1:P.qqSample:end,:)));
        mrAmp = bsxfun(@rdivide, mrAmp, std(mrAmp,0,2));
        Sevt.trFet = shiftdim(mrAmp, -1);            
        if P.vcFet(end)=='2', Sevt.trFet=Sevt.trFet.^2; end
    case 'wavcov'
        Sevt.trFet = shiftdim(getWavCov(mrWav, Sevt, P), -1);
    case {'slope', 'slope2'}
        Sevt.trFet = shiftdim(getSlope_(mrWav, Sevt.viSpk, P), -1);
        if P.vcFet(end)=='2', Sevt.trFet=Sevt.trFet.^2; end    
%     case 'slope2' %slope squared
%         Sevt.trFet = shiftdim(getSlope(mrWav, Sevt.viSpk, P).^2, -1);
%         Sevt.trFet = getSlope2(mrWav, Sevt.viSpk, P);
    case 'area' %area under the curve
        Sevt.trFet = shiftdim(getArea(mrWav, Sevt.viSpk, P), -1);
    case 'diff248'
        Sevt.trFet = fet_diff248_(mrWav, Sevt.viSpk, P);
    otherwise
        error(['invalid feature: ', P.vcFet]);
end
Sevt.dimm_fet = size(Sevt.trFet);
end %func


function Sevt = load_evt_(vcFile_prm)
Sevt = load_evt(vcFile_prm);
end %func


function save_evt_(vcFile_evt, Sevt)
% save Sevt info
S0 = get(0, 'UserData');
S0.Sevt = Sevt;
set(0, 'UserData', S0);

if iscell(Sevt.trFet)
    save_struct_(vcFile_evt, Sevt);
    return;
end
trFet = Sevt.trFet;
Sevt.dimm_fet = size(trFet);
Sevt.trFet = [];
save_struct_(vcFile_evt, Sevt);

% save fet file as a binary file
if ~isempty(trFet)
    vcFile_fet = strrep(vcFile_evt, '_evt.mat', '.fet');
    write_bin_(vcFile_fet, trFet);    
%     fprintf('Saving to %s...', vcFile_fet);   t1=tic;
%     fid = fopen(vcFile_fet, 'W');
%     fwrite(fid, trFet, class(trFet));
%     fclose(fid);
%     fprintf(' took %0.1fs\n', toc(t1));
end
end %func


function mrCov = getWavCov(mrWav, Sevt, P)
fMaxRef = 1;
fMeanRef = 0;

miSites = findNearSites(P.mrSiteXY, P.maxSite, P.viSiteZero);
nSpk = numel(Sevt.viSpk);
nSites = size(mrWav,2);
mrCov = zeros(nSites, nSpk, 'single');
for iSpk=1:nSpk
    iTime1 = Sevt.viSpk(iSpk);    
    iSite1 = Sevt.viSite(iSpk);
%     maxAmp1 = single(mrWav(iTime1, iSite1));
    mrWav1 = squeeze(mr2tr2(mrWav, iTime1, P.spkLim));    
    mrWav1 = single(mrWav1);
%     mrWav1 = bsxfun(@minus, mrWav1, mean(mrWav1));
    if fMaxRef
        vrRef = mrWav1(:,iSite1);        
        mrWav1 = bsxfun(@minus, mrWav1, vrRef);
    elseif fMeanRef
        vrRef = mean(mrWav1(:, miSites(:,iSite1)), 2);        
        mrWav1 = bsxfun(@minus, mrWav1, vrRef);
    end    
    
    mrCov(:,iSpk) = sqrt(sum(mrWav1.^2));
%     mrCov(:,iSpk) = sum(mrWav1);
%     mrWav1 = bsxfun(@times, mrWav1, mrWav1(:,iSite1));
%     mrCov(:,iSpk) = sum(mrWav1);
end
end %func


function Sevt = loadEvtCache_(vcFile_evt)
% Uses Sevt from workspace
try
    Sevt = evalin('base', 'Sevt');
    vcFile_evt_cache = subsFileExt(Sevt.P.vcFile_prm, '_evt.mat');
    if ~strcmpi(vcFile_evt_cache, vcFile_evt)
        Sevt = []; %different file
        return;
    else
        fprintf('Loaded %s from cache\n', vcFile_evt);
    end
catch
    Sevt = [];
%     disp(lasterror());
    return;
end
end %func


function mrEnergy = getEnergy_(mrWav, viSpk, P)

fprintf('Calculating energy...'); t1=tic;
nChans = size(mrWav, 2);
mrEnergy = zeros(numel(viSpk), nChans, 'single');
for iCh=1:nChans
    mrWav1 = vr2mr2(single(mrWav(:,iCh)), viSpk, P.spkLim);
    mrEnergy(:,iCh) = std(mrWav1, 1);
end
mrEnergy = mrEnergy' * P.uV_per_bit;
fprintf(' took %0.1fs\n', toc(t1));  
end


function mrEnergy = getEnergy_fft_(mrWav, viSpk, P)

fprintf('Finding energy density...'); t1=tic;

mrWav1 = squeeze(mr2tr2(mrWav, viSpk, P.spkLim));    
trWav1_dim = size(mrWav1);
mrWav1 = reshape(mrWav1, trWav1_dim(1), []); %convert to 2d
mrWav1 = single(mrWav1);
mrWav1 = bsxfun(@minus, mrWav1, mean(mrWav1));
mrWav1 = bsxfun(@times, mrWav1, hamming(trWav1_dim(1))); %apply hamming window
mrWav1 = abs(fft(mrWav1));
mrWav1 = mrWav1(1:ceil(end/2),:);
mrEnergy = squeeze(mean(mrWav1.^2)); %only use half for speed
mrEnergy = reshape(mrEnergy, trWav1_dim(2), []);

fprintf(' took %0.1fs\n', toc(t1));

% if fSubtDc
%     trWav1 = trWav1(2:ceil(end/2),:,:); %no dc power
% else
%     trWav1 = trWav1(1:ceil(end/2),:,:); %no dc power
%     trWav1(1,:,:) = trWav1(1,:,:)/2; %half
% end

% subtract high-freq noise
% if fSubtHighFreq
%     trWav1_dim = size(mrWav1);
%     mrWav1 = reshape(mrWav1, trWav1_dim(1), []);
%     mrWav1 = bsxfun(@minus, mrWav1(1:end-1,:), mrWav1(end,:));
%     trWav1_dim(1) = trWav1_dim(1)-1;
%     mrWav1 = reshape(mrWav1, trWav1_dim);
% end
% if fCutHalf
%     mrWav1 = mrWav1(1:round(end/2),:,:);
% end

    
end


function mrSlope = getSlope_(mrWav, viSpk, P)
fprintf('Finding slope...'); t1=tic;

nMax = size(mrWav,1);
viSpk1 = saturateLim(viSpk + P.slopeLim(1), 1, nMax);
viSpk2 = saturateLim(viSpk + P.slopeLim(end), 1, nMax);
mrSlope = (mrWav(viSpk1, :) - mrWav(viSpk2, :))';
mrSlope = single(mrSlope) * P.uV_per_bit; %nSites x nSpk

fprintf(' took %0.1fs\n', toc(t1));
end


function mrArea = getArea(mrWav, viSpk, P)
if nargin<3, P.uV_per_bit=1; end
fprintf('Finding area...'); t1=tic;

viSpk1 = viSpk - P.nSlope;
viSpk1(viSpk1<1) = 1;
mrArea = single(mrWav(viSpk1, :)) * (-P.nSlope);
nWav = size(mrWav,1);
for i = 1:P.nSlope
    viSpk1 = viSpk1 + 1;    
    viSpk1(viSpk1>nWav)=nWav;
    mrArea = mrArea + single(mrWav(viSpk1, :));
end
mrArea = mrArea' * P.uV_per_bit;

fprintf(' took %0.1fs\n', toc(t1));
end


function trSlope = getSlope2(mrWav, viSpk, P)
% two slopes, rising and falling
if nargin<3, P.uV_per_bit=1; end

fprintf('Finding dual slope...'); t1=tic;
nWav = size(mrWav,1);
viSpk1 = viSpk - P.nSlope;  viSpk1(viSpk1<1) = 1;
viSpk2 = viSpk - round(P.nSlope/2);  viSpk2(viSpk2<1) = 1;
mrWav0 = mrWav(viSpk,:)';
mrWav1 = mrWav(viSpk1,:)';
mrWav2 = mrWav(viSpk2,:)';
trSlope = cat(3, mrWav0 - mrWav1, mrWav2 - mrWav1);
trSlope = permute(trSlope, [3 1 2]);
trSlope = single(trSlope) * P.uV_per_bit; %in uV unit
fprintf(' took %0.1fs\n', toc(t1));
end


% function process_mutlfile(fid, P)
% Merge multiple files and process as single
% %     vcDir1 = sprintf('%s%s\\', vcDir, csAnimals{iAnimal});
% 
% error('not implemented yet. cluster multiple data files'); % merge data files and cluster together
% 
% csFiles = P.vcFile; %backup
% 
% %recursive execution
% for iFile=1:numel(csFiles)
%     vcFile = csFiles{iFile};
%     P.vcFile = vcFile;
%     ffprintf([1,fid], vcFile);        
%     if ~P.fRun, continue; end %just list files to be processed        
%     [vcFile_spk, vcFile_evt, vcFile_clu] = subsFileExt(vcFile, '.spk', '_evt.mat', '_clu.mat');
% 
%     if ~exist(vcFile_clu, 'file') || P.fOverwriteClu
%         % event detection
%         if ~exist(vcFile_evt, 'file') || P.fOverwriteEvt
%             fRunEvt = 1;
%         else
%             try
%                 load(vcFile_evt);
%                 fRunEvt = 0;
%             catch
%                 fRunEvt = 1;
%             end
%         end
%         if fRunEvt
%             if ~exist(vcFile_spk, 'file') || P.fOverwriteSpk
%                 % load file. not saving it
%                 mrWav = export_spk(vcFile, 'freqLim', P.freqLim, 'fSave', P.fSaveSpk);
%                 P.nChans = size(mrWav,2);
%                 if P.fSaveSpk, ffprintf([1,fid], ['wrote to ', vcFile_spk]); end
%             else
%                 mrWav = loadWavFile(vcFile_spk, P.nChans);
%             end
%             Sevt = jrclust_event(mrWav, P);
%             clear mrWav; %don't need anymore
%             writeToMat(vcFile_evt, Sevt);
%             ffprintf([1,fid], ['wrote to ', vcFile_evt]); 
%         end            
%     end
% end %for file
% %run cluster
% if P.fRunCluster
%     Sclu = jrclust_sort(Sevt, P);
%     clear Sevt; %don't need after
%     writeToMat(vcFile_clu, Sclu);
%     ffprintf([1,fid], ['wrote to ', vcFile_clu]);
% end
% end %func


% function [viSite2Chan, mrSiteXY, vrSiteHW] = read_prb_file_(vcFilePrb)
% P_ = file2struct(vcFilePrb);
% viSite2Chan = P_.channels;
% mrSiteXY = P_.geometry;
% vrSiteHW = P_.pad;
% end

function P = read_meta_file_(vcFile_meta)
P = [];
if ~exist(vcFile_meta, 'file'), return; end
try
    Smeta = read_whisper_meta(vcFile_meta);
    if isempty(Smeta), return; end
    P = struct('sRateHz', Smeta.sRateHz, 'uV_per_bit', Smeta.scale, 'nChans', Smeta.nChans, 'probe_file', [Smeta.vcProbe, '.prb'], 'vcDataType', Smeta.vcDataType);
    P.Smeta = Smeta;
catch
    disp(lasterror());
end
end %func


function [P, vcPrompt] = create_prm_file_(vcFile_bin, vcFile_prb, vcFile_template, fAsk)
if nargin<2, vcFile_prb = ''; end
if nargin<3, vcFile_template = ''; end
if nargin<4, fAsk = 1; end
P0 = file2struct('default.prm');  %P = defaultParam();
if any(vcFile_bin=='*') %merge multiple files
    vcFile_bin1 = vcFile_bin;
    vcFile_bin = strrep(vcFile_bin1, '*', 'all');
    merge_binfile_(vcFile_bin, vcFile_bin1);
end
if ~exist(vcFile_bin, 'file')
    P = []; 
    vcPrompt = sprintf('%s does not exist.\n', vcFile_bin);    
    fprintf(2, '%s\n', vcPrompt); return;
end
if matchFileExt(vcFile_template, '.prm') %template file provided
    P.template_file = vcFile_template;
end
% append probe file
if ~isempty(vcFile_prb)
    vcPostfix = ['_', subsFileExt(vcFile_prb, '')];
    P.vcFile_prm = subsFileExt(vcFile_bin, [vcPostfix, '.prm']);
    P.probe_file = vcFile_prb;
else
    P.vcFile_prm = subsFileExt(vcFile_bin, '.prm');
end
if exist(P.vcFile_prm, 'file') && fAsk
    vcAns = questdlg('File already exists. Overwrite prm file?', 'Warning', 'Yes', 'No', 'No');
    if ~strcmpi(vcAns, 'Yes')
        P = [];
        vcPrompt = 'Cancelled by user.';
        return;
    end
end

% Load meta file
[~,~,vcExt] = fileparts(vcFile_bin);
switch lower(vcExt)
    case '.ncs' %neuralynx format
        [~, P_meta, P.vcFile] = ncs2dat(vcFile_bin, 'viChan_bin', P.viChan_bin, 'tlim_load', P.tlim_load);
    case {'.bin', '.dat'}
        P.vcFile = vcFile_bin;
        P_meta = read_meta_file_(subsFileExt(P.vcFile, '.meta'));
        if ~isempty(vcFile_prb) && ~isempty(P_meta)
            P_meta = rmfield(P_meta, 'probe_file');
        end
end

% Load prb file
if isfield(P, 'template_file')
    P = appendStruct(file2struct(P.template_file), P);
end

% if ~isempty(P_meta)
%     disp('Meta file does not exist.'); 
%     P.vcDataType = P0.vcDataType;  
%     S_prb = file2struct(P.probe_file);
%     P.nChans = numel(S_prb.channels);
%     P.sRateHz = P0.sRateHz;
% else
P = appendStruct(P, P_meta);    
P = appendStruct(P, file_info_(vcFile_bin));

P.duration_file = P.nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans / P.sRateHz; %assuming int16

% if ~isempty(vcProbeType), P.probe_file = vcProbeType; end
copyfile('default.prm', P.vcFile_prm, 'f');
edit_prm_file_(P, P.vcFile_prm);
vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
disp(vcPrompt);
if fAsk, edit(P.vcFile_prm); end % Show settings file
end %func


function S1 = copyStruct_(S1, S2, csFields)
for i=1:numel(csFields)
    if isfield(S2, csFields{i}), S1.(csFields{i}) = S2.(csFields{i}); end
end
end %func


function S = remove_struct_(S, varargin)
% remove fields from a struct
for i=1:numel(varargin)
    if isfield(S, varargin{i})
        S = rmfield(S, varargin{i});
    end
end
end %func


function P = load_prb_(vcFile_prb, P)
% append probe file to P
P.probe_file = vcFile_prb;
%     [P.viSite2Chan, P.mrSiteXY, P.vrSiteHW, P.cviShank] = read_prb_file(vcFile_prb);
S_prb = file2struct(vcFile_prb);
P.viSite2Chan = S_prb.channels;
P.mrSiteXY = S_prb.geometry;
P.vrSiteHW = S_prb.pad;
S_prb = remove_struct_(S_prb, 'channels', 'geometry', 'pad', 'ref_sites', 'viHalf', 'i', 'vcFile_file2struct');
if isfield(S_prb, 'viShank')
    P.viShank = S_prb.viShank;
else
    P.viShank = ones(1, numel(P.viSite2Chan));
end
% P = copyStruct_(P, S_prb, {'cviShank', 'maxSite', 'um_per_pix'});
P.viChan_aux = setdiff(1:P.nChans, 1:max(P.viSite2Chan)); %aux channel. change for
P = appendStruct(P, S_prb);
end %func


function P = file_info_(vcFile)
S_dir = dir(vcFile);
P.vcDate_file = S_dir.date;
P.nBytes_file = S_dir.bytes;
end


function [P, vcFile_prm] = loadParam_(vcFile_prm)
% Load prm file
P0 = file2struct('default.prm');  %P = defaultParam();
%P0 = appendStruct(P0, file2struct('./template.prm')); %overwrite template
P = file2struct(vcFile_prm);
if ~isfield(P, 'template_file'), P.template_file = ''; end
if ~isempty(P.template_file)
    P0 = appendStruct(P0, file2struct(P.template_file));
end
P.vcFile_prm = vcFile_prm;
% todo: substitute bin file path
if ~exist(P.vcFile, 'file')
    P.vcFile = replacePath(P.vcFile, vcFile_prm);
    if ~exist(P.vcFile, 'file'), fprintf(2,'file does not exist.\n'); return; end
end
% Load prb file
if ~isfield(P, 'probe_file'), P.probe_file = P0.probe_file; end
try    
    if ~exist(P.probe_file, 'file')
        P.probe_file = replacePath(P.probe_file, vcFile_prm); 
        if ~exist(P.probe_file, 'file'), error('prb file does not exist'); end
    end
    P0 = load_prb_(P.probe_file, P0);
catch
    fprintf('loadParam: %s not found.\n', P.probe_file);
end

% Load prb file
P = appendStruct(P0, P);    

% computed fields
if isempty(P.vcFile_prm), P.vcFile_prm = subsFileExt(P.vcFile, '.prm'); end
P.spkRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);
P.spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
if isempty(P.spkLim_ms_fet), P.spkLim_ms_fet = P.spkLim_ms; end
P.spkLim_fet = round(P.spkLim_ms_fet * P.sRateHz / 1000);
P.slopeLim = round(P.slopeLim_ms * P.sRateHz / 1000);
if ~isempty(P.nDiff_ms_filt)
    P.nDiff_filt = ceil(P.nDiff_ms_filt * P.sRateHz / 1000);
else
    P.nDiff_filt = 0;
end
try P.miSites = findNearSites(P.mrSiteXY, P.maxSite, P.viSiteZero); catch; end %find closest sites
P.sRateHz_lfp = P.sRateHz / P.nSkip_lfp;        %LFP sampling rate
P.bytesPerSample = bytesPerSample_(P.vcDataType);
if isempty(P.vcFile_prm), P.vcFile_prm = subsFileExt(P.vcFile, '.prm'); end %backward compatibility
if isfield(P, 'gain_boost'), P.uV_per_bit = P.uV_per_bit / P.gain_boost; end
try 
    if isempty(P.maxSite_fet), P.maxSite_fet = P.maxSite; end
    if isempty(P.maxSite_detect), P.maxSite_detect = P.maxSite; end
    if isempty(P.maxSite_sort), P.maxSite_sort = P.maxSite; end
    if isempty(P.maxSite_pix), P.maxSite_pix = P.maxSite; end
    if isempty(P.maxSite_dip), P.maxSite_dip = P.maxSite; end
    if isempty(P.maxSite_merge), P.maxSite_merge = P.maxSite; end
catch
    disp(lasterr);
end
%set0(P); %put in the memory
edit(P.vcFile_prm); % Show settings file
% delete('temp_jrclust_*.m'); %in case not deleted
end


function n = bytesPerSample_(vcDataType)
switch lower(vcDataType)
    case {'int16', 'uint16'}
        n = 2;
    case {'single', 'float', 'int32', 'uint32'}
        n = 4;
    case {'double', 'int64', 'uint64'}
        n = 8;
end
end %func


function cs = first_string_(cs)
% extract first string
for i=1:numel(cs)
    if isempty(cs{i}), continue; end
    cs1 = textscan(cs{i}, '%s', 'Delimiter', {' ','='});
    cs1 = cs1{1};
    cs{i} = cs1{1};
end
end %func


function edit_prm_file_(P, vcFile_prm)
% turn a struct to file
csLines = file2cellstr(vcFile_prm); %read to cell string
csLines_var = first_string_(csLines);

csName = fieldnames(P);
csValue = cellfun(@(vcField)P.(vcField), csName, 'UniformOutput',0);
for i=1:numel(csName)
    vcName = csName{i}; %find field name with 
    if isstruct(csValue{i}), continue; end %do not write struct
    vcValue = field2str(csValue{i});
    iLine = find(strcmpi(csLines_var, vcName));
    if numel(iLine)>1
        error(['edit_prm_file_ error: ' vcName]); 
    elseif isempty(iLine) %did not find, append
        csLines{end+1} = sprintf('%s = %s;', vcName, vcValue);
    else    
        vcComment = getCommentExpr_(csLines{iLine});        
        csLines{iLine} = sprintf('%s = %s;\t\t\t%s', vcName, vcValue, vcComment);    
    end
end
cellstr2file(vcFile_prm, csLines);
end %func


function vcComment = getCommentExpr_(vcExpr)
iStart = strfind(vcExpr, '%');
if isempty(iStart), vcComment = ''; return; end
vcComment = vcExpr(iStart(1):end);
end %func


% function write_prm_file_(P, vcFile_prm)
% % edit a parameter file
% 
% csFieldname = fieldnames(P);
% 
% fid = fopen(vcFile_prm, 'w');
% for i=1:numel(csFieldname)
%     val = getfield(P, csFieldname{i});
%     fprintf(fid, '%s = %s;\n', csFieldname{i}, field2str(val));
% end %for
% fclose(fid);
% fprintf('Wrote to %s\n', vcFile_prm);
% end %func


function [csFiles, vlDir] = getFileList_(vcFile0, vcExt)
% returns a file list containing the extension provided, and whether it's
% ext
[vcDir0, ~, vcExt0] = fileparts(vcFile0);
if any(strcmpi(vcExt0, vcExt))
    csFiles = {vcFile0};
    vlDir = 0;
    return;
end
dir0 = dir(vcFile0);
csFiles = {dir0.name};
vlDir = cell2mat({dir0.isdir});

% remove mismatches
vlKill = ismember(csFiles, {'.', '..'}) | ~matchFileExt(csFiles, vcExt, vlDir);
csFiles(vlKill) = [];
vlDir(vlKill) = [];

% make full path
csFiles = cellfun(@(x)[vcDir0, filesep(), x], csFiles, 'UniformOutput', 0);

% sort by creation time

end %func


function plotLfp_jrclust_(vcFile_prm)
% @TODO
if ~matchFileExt(vcFile_prm, '.prm')
    vcFile_prm = subsFileExt(vcFile_prm, '.prm');
end
P = loadParam_(vcFile_prm);

vcFile_lfp = subsFileExt(vcFile_prm, '.lfp');
P1.maxAmp = P.maxAmp_lfp;
P1.sRateHz = P.sRateHz_lfp;
P1.freqLim = [];
P1.tlim = P.tlim_lfp;
% P1.fFilt = 0;
P1.nChans = numel(P.viSite2Chan);
P1.nSites_ref = P.nSites_ref;
jrclust_traces(vcFile_lfp, P1);

end %func


function plotRef_jrclust_(vcFile_prm)
if ~matchFileExt(vcFile_prm, '.prm')
    vcFile_prm = subsFileExt(vcFile_prm, '.prm');
end
P = loadParam_(vcFile_prm);

fid = fopen(subsFileExt(P.vcFile, '.ref'), 'r');    
vrWav_ref = fread(fid, inf, '*int16');
fclose(fid);
figure; plot(vrWav_ref);
set(gcf, 'Name', vcFile_prm);
end %func


function plotTraces_jrclust_(vcFile_prm)
if ~matchFileExt(vcFile_prm, '.prm')
    fprintf(2, 'Must provide .prm file\n');
    return;
end
% S0 = get(0, 'UserData');
% % try P = S0.P; catch; end
if ~exist(vcFile_prm, 'file'), fprintf(2, '.prm file does not exist.'); return; end
P = loadParam_(vcFile_prm);
Sclu = load_clu_(vcFile_prm);
Sevt = load_evt_(vcFile_prm);
% S0 = appendStruct(S0, makeStruct(P, Sclu, Sevt));
% set(0, 'UserData', S0);

if ~isempty(Sevt)
    P.viSpk = Sevt.viSpk;     
    P.viSite = Sevt.viSite;
end
if ~isempty(Sclu), P.viClu = Sclu.viClu; end
if isempty(Sevt)
    try
        S_gt = load(subsFileExt(vcFile_prm, '_gt.mat'));
        uiwait(msgbox('Spikes loaded from ground truth. Press OK to continue.'));
        P.viSpk = S_gt.viTime;
        P.viSite = ones(size(P.viSpk));
        P.viClu = S_gt.viClu;
    catch
        disp('No ground truth files');
    end
end

% pass .bin file
vcFile_bin = P.vcFile;
if ~exist(vcFile_bin, 'file')
    vcFile_bin = replacePath(vcFile_bin, vcFile_prm);
    if ~exist(vcFile_bin, 'file')
        fprintf(2,'%s not found\n', vcFile_bin);
        return;
    end
end

jrclust_traces(vcFile_bin, P);
end %func


function vcFile = filename_append(vcFile, vcEnd)
[vcDir, vcFname, vcExt] = fileparts(vcFile);
if isempty(vcDir), vcDir='.'; end
vcFile = [vcDir, filesep(), vcFname, vcEnd, vcExt];
end %func


function vcFile = filename_format(vcFile)
% add proper irectory and extensions if missing
[vcDir, vcFname, vcExt] = fileparts(vcFile);
if isempty(vcDir), vcDir='.'; end
if isempty(vcExt) && numel(vcFname)>3
    switch lower(vcFname(end-3:end))
        case {'_clu', '_evt'}
            vcExt = '.mat';
        case {'_prm', '_prb'}
            vcExt = '.m';    
    end
else
    vcExt = '';
end
vcFile = [vcDir, filesep(), vcFname, vcExt];
end %func


function vrWav1 = spatial_subtract_(mrWav, iChan, P)
viSites1 = P.miSites(2:P.maxSite*2+1, iChan);
vrWav1 = mrWav(:,iChan) - int16(mean(mrWav(:,viSites1),2));
end


function gvrSum = spatial_average_gpu_(mnWav, iChan, P)
% GPU cache building
persistent cgvrWav
if nargin==0, cgvrWav=[]; return; end
if isempty(P.maxSite_detect), P.maxSite_detect = P.maxSite; end
nChan_spatial = P.maxSite_detect * 2 + 1;
nChans = size(mnWav,2);
if isempty(cgvrWav), cgvrWav = cell(nChans, 1); end
% gvrSum = zeros(size(mnWav,1), 1, 'single', 'gpuArray');
gvrSum = zeros(size(mnWav,1), 1, 'int16', 'gpuArray');
viSites1 = P.miSites(1:nChan_spatial, iChan); 
vlChan_sum = ismember(1:nChans, viSites1);
for iChan1 = 1:nChans
    if vlChan_sum(iChan1)
        if isempty(cgvrWav{iChan1}) % selective gpu upload
            cgvrWav{iChan1} = gpuArray(mnWav(:, iChan1));        
        end
%         gvrSum = gvrSum + single(cgvrWav{iChan1}); % optional divide
            gvrSum = gvrSum + cgvrWav{iChan1}; % optional divide
    else
        cgvrWav{iChan1} = []; %free from memory
    end    
end
gvrSum = single(gvrSum) / numel(viSites1);
% gvrSum = int16(gvrSum);
end %func


function vrWav1 = spatial_average_(mnWav, iChan, P)
nChan_spatial = P.maxSite_detect * 2 + 1;
viSites1 = P.miSites(1:nChan_spatial, iChan); 
% vrWav1 = zeros(size(mnWav,1), 1, 'single');
% tic
for i=1:numel(viSites1)
    if i==1
        vrWav1 = mnWav(:, viSites1(i));
    else
        vrWav1 = vrWav1 + mnWav(:, viSites1(i));
    end
end
if isa(mnWav, 'int16')
    vrWav1 = int16(single(vrWav1) / numel(viSites1));
elseif isa(mnWav, 'double') || isa(mnWav, 'single')
    vrWav1 = single(vrWav1) / numel(viSites1);
elseif isa(mnWav, 'int32')
    vrWav1 = int32(single(vrWav1) / numel(viSites1));
else
    error('spatial_average_: unsupported class');
end
end %func


function gvrSum = spatial_csd_gpu_(mnWav, iChan, P)
% GPU cache building. returns negative value for negative peak

persistent cgvrWav
if nargin==0, cgvrWav=[]; return; end
nChan_spatial = P.maxSite * 2 + 1;
nChans = size(mnWav,2);
if isempty(cgvrWav), cgvrWav = cell(nChans, 1); end
gvrSum = zeros(size(mnWav,1), 1, 'single', 'gpuArray');
viSites1 = P.miSites(1:nChan_spatial, iChan); 
vlChan_sum = ismember(1:nChans, viSites1);
if isempty(cgvrWav{iChan}), cgvrWav{iChan} = gpuArray(mnWav(:, iChan)); end
for iChan1 = 1:nChans
    if vlChan_sum(iChan1)
        if isempty(cgvrWav{iChan1}) % selective gpu upload
            cgvrWav{iChan1} = gpuArray(mnWav(:, iChan1));        
        end
        gvrSum = gvrSum + single(abs(cgvrWav{iChan} - cgvrWav{iChan1})); % optional divide
    else
        cgvrWav{iChan1} = []; %free from memory
    end    
end
% gvrSum = rms_filt(gvrSum, P.nRms_filt) * (-1 / numel(viSites1));
% gvrSum = gvrSum / numel(viSites1);
% gvrSum = single(gvrSum) / numel(viSites1);
end %func


function vrWav1 = spatial_imec2_(mrWav, iChan, P)
viSites1 = P.miSites(1:4, iChan); %use two immediate neighbors
vrMean = single([4 2 1 1]' / 8); %matrix product
vrWav1 = single(mrWav(:,viSites1)) * vrMean;
end


function [cviSpk, cvrSpk, vrThresh, cmrSpk] = spikeDetect_(mrWav, P)
% P: spkLim = [-10 24], qqFactor = 4.5, P.qqSample = 4, P.minAmpLimit = []
if ~isfield(P, 'vrPv'), P.vrPv = []; end
if ~isfield(P, 'nlim_detect'), P.nlim_detect = []; end

fprintf('Detecting spikes\n'); t1=tic;
[nSamples, nChan] = size(mrWav);
[cviSpk, cvrSpk] = deal(cell(nChan,1));
vrThresh = zeros(nChan, 1);
if nargout>=4
    fSpkTable = 1;
    cmrSpk = cell(nChan,1); 
else
    fSpkTable = 0;
end

if isempty(P.viSiteZero)
    P.viSiteZero = find(max(mrWav) == min(mrWav));
end
viChan = setdiff(1:nChan, P.viSiteZero);
spatial_filter_(); %clear cache
% P.fGpu = 1; %parfor used so do this in main memory
% nRms_filt = round(P.rms_filt_ms * P.sRateHz / 1000);
% if ~P.fRms_detect, nRms_filt = 0; end
% if P.nDiff_filt > 0
%     vi0 = gpuArray(int32(1:nSamples));
%     cvi0 = {min(max(vi0+P.nDiff_filt+1, 1), nSamples), ...
%             min(max(vi0-P.nDiff_filt, 1), nSamples), ...
%             min(max(vi0+P.nDiff_filt, 1), nSamples), ...
%             min(max(vi0-P.nDiff_filt-1, 1), nSamples)};
% end
% P.fGpu=1;  
if isempty(P.nDiff_filt), P.nDiff_filt = 0; end
for iCh = viChan
    try
        if isempty(P.nDiff_filt), P.nDiff_filt=0; end    
    %     vrWav1 = matched_filt_(spatial_filter_(mrWav, iCh, P), P.vrPv, P.spkLim);    
        if P.fGpu
            vrWav1=gpuArray(mrWav(:,iCh)); 
        else
            vrWav1 = mrWav(:,iCh);
        end
    %     if nRms_filt>0, vrWav1 = -1 * rms_filt(vrWav1, nRms_filt); end
        if ~isa(vrWav1,'int16'), vrWav1 = int16(vrWav1); end
        [cviSpk{iCh}, cvrSpk{iCh}, vrThresh(iCh)] = spikeDetectSingle_fast_(vrWav1, P);
        if ~isempty(P.vrPv)
            cviSpk{iCh} = cviSpk{iCh} - P.spkLim(2); 
        end %matched filter
        fprintf('.');
    catch
        disp(lasterr)
    end
end

% refractory period
parfor iCh=1:nChan
    [cviSpk{iCh}, cvrSpk{iCh}] = spike_refrac(cviSpk{iCh}, cvrSpk{iCh}, P.spkRefrac);
end
spatial_filter_(); %clear cache
nSpks = sum(cellfun(@numel, cviSpk));
fprintf('\n\tDetecting %d spikes took %0.1fs.\n', nSpks, toc(t1));
% vrThresh = cell2mat_(vrThresh);
end %func


function vrWav1 = spatial_filter_(mrWav, iChan, P)
% reset gpu cache
if nargin==0, spatial_csd_gpu_(); spatial_average_gpu_(); return; end

switch lower(P.vcSpatialFilter) %set to 'none' to skip
    case 'csd'
        vrWav1 = spatial_csd_gpu_(mrWav, iChan, P); %check for GPU use
    case 'subtract'
        vrWav1 = spatial_subtract_(mrWav, iChan, P);
    case {'average', 'mean'}
        if P.fGpu
            vrWav1 = spatial_average_gpu_(mrWav, iChan, P); %check for GPU use
        else
            vrWav1 = spatial_average_(mrWav, iChan, P); %check for GPU use
        end
    case 'imec2'
        vrWav1 = spatial_imec2_(mrWav, iChan, P);
    otherwise %spatial filter
        if isempty(P.nlim_detect)
            vrWav1 = mrWav(:,iChan); 
        else
            vrWav1 = mrWav(P.nlim_detect(1):P.nlim_detect(end), iChan); % time range specified
        end
        if P.fGpu
%             vrWav1 = single(gpuArray(vrWav1));
            vrWav1 = (gpuArray(vrWav1));
        else
            vrWav1 = single(vrWav1);
        end
end 
end %func


function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_fast_(vrWav1, P)
% P: spkThresh, qqSample, qqFactor, fGpu, uV_per_bit
% vrWav1 can be either single or int16

% Determine threshold
MAX_SAMPLE_QQ = 300000; 
if isempty(P.spkThresh)
    if P.nDiff_filt>0
        thresh1 = median(abs(subsample_vr(vrWav1, MAX_SAMPLE_QQ)));
        thresh1 = int16(single(thresh1)* P.qqFactor / 0.6745);
    else %uncentered
        vrWav1_sub = subsample_vr(vrWav1, MAX_SAMPLE_QQ);
        med1 = median(vrWav1_sub);
        thresh1 = median(abs(vrWav1_sub - med1));
        thresh1 = int16(single(thresh1)* P.qqFactor / 0.6745) + abs(med1);
    end
else
    thresh1 = int16(P.spkThresh); 
end

% detect valley turning point. cannot detect bipolar
% pick spikes crossing at least three samples
viSpk1 = find_peak_(vrWav1, thresh1);
if P.fDetectBipolar
   viSpk1 = [viSpk1; find_peak_(-vrWav1, thresh1)]; 
   viSpk1 = sort(viSpk1);
end
if isempty(viSpk1)
    viSpk1 = double([]);
    vrSpk1 = int16([]);
else
    vrSpk1 = vrWav1(viSpk1);
    % Remove spikes too large
    if ~isempty(P.spkThresh_max_uV)
        thresh_max1 = int16(abs(P.spkThresh_max_uV) / P.uV_per_bit);
    %     viA1 = find(vrSpk1 > -thresh_max1);
        viA1 = find(abs(vrSpk1) < abs(thresh_max1));
        viSpk1 = viSpk1(viA1);
        vrSpk1 = vrSpk1(viA1); 
    end    
end

if P.fGpu
    [viSpk1, vrSpk1, thresh1] = multifun(@gather, viSpk1, vrSpk1, thresh1);
end
% thresh1 = int16(thresh1);
end %func


function viSpk1 = find_peak_(vrWav1, thresh1)
viSpk1 = [];
vi2 = find(vrWav1 < -thresh1);
if isempty(vi2), thresh1 = 0; return; end
viSpk1 = vi2(find(diff(diff(vrWav1(vi2))>0)>0) + 1); % only negative peak
if isempty(viSpk1), thresh1 = 0; return; end
if viSpk1(1) <= 1, viSpk1(1) = 2; end
if viSpk1(end) >= numel(vrWav1), viSpk1(end) = numel(vrWav1)-1; end
viSpk1 = viSpk1(vrWav1(viSpk1-1) < -thresh1 & vrWav1(viSpk1+1) < -thresh1);
end %func


function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_(vrWav1, P)
% P: spkThresh, qqSample, qqFactor, fGpu, uV_per_bit
% vrWav1 can be either single or int16
MAX_SAMPLE_QQ = 300000; 

if P.fGpu
    try vrWav1 = gpuArray(vrWav1); catch; end
end
if isempty(P.spkThresh)
    thresh1 = median(abs(subsample_vr(vrWav1, MAX_SAMPLE_QQ)));
%     thresh1 = median(abs(vrWav1(1:P.qqSample:end)));
    thresh1 = int16(single(thresh1)* P.qqFactor / 0.6745);
else
    thresh1 = int16(P.spkThresh); 
end

% detect valley turning point
if ~P.fDetectBipolar
    viSpk1 = find(diff(diff(vrWav1)>0)>0) + 1; % only negative peak
    vrSpk1 = vrWav1(viSpk1);
else
    viSpk1 = find(diff(diff(vrWav1)>0)) + 1; % negative and positive
    vrSpk1 = vrWav1(viSpk1);
    vlPos = vrSpk1 > 0;  
    vrSpk1(vlPos) = -vrSpk1(vlPos); %make positive negative
end

% use below threshold spikes only
if isempty(P.spkThresh_max_uV)
    viA1 = find(vrSpk1 < -thresh1);
else    
    thresh_max1 = int16(abs(P.spkThresh_max_uV) / P.uV_per_bit);
    viA1 = find(vrSpk1 < -thresh1 & vrSpk1 > -thresh_max1);
end
viSpk1 = viSpk1(viA1);
vrSpk1 = vrSpk1(viA1); 

if P.fGpu
    viSpk1 = gather(viSpk1);
    vrSpk1 = gather(vrSpk1);
    thresh1 = gather(thresh1);
end

% refractory period
[viSpk1, vrSpk1] = spike_refrac(viSpk1, vrSpk1, P.spkRefrac*2); %its own chan has 2x refrac
if isempty(viSpk1)
    viSpk1 = double([]);
    vrSpk1 = single([]);
end
end %func


function [viSpk, vrSpk, viSite] = spikeMerge_(cviSpk, cvrSpk, P)
% provide spike index (cviSpk) and amplitudes (cvrSPk) per sites
if nargin < 3
    P = struct('spkLim', [-10 24], 'maxSite', 2.5);    
end

nSites = numel(cviSpk);
viSpk = cell2mat_(cviSpk);      vrSpk = cell2mat_(cvrSpk);
viSite = cell2mat_(cellfun(@(vi,i)repmat(i,size(vi)), cviSpk, num2cell((1:nSites)'), 'UniformOutput', false));
[viSpk, viSrt] = sort(viSpk);   vrSpk = vrSpk(viSrt);   viSite = viSite(viSrt);
viSite = int32(viSite); 
viSpk = int32(viSpk);   

% jitter = ceil(diff(spkLim)/2);
% remove multiple detection from other channels
[cviSpkA, cvrSpkA, cviSiteA] = deal(cell(nSites,1));
% if isempty(P.maxSite_detect)
%     P.miSites_detect = P.miSites;
% else
%     P.miSites_detect = findNearSites(P.mrSiteXY, P.maxSite_detect, P.viSiteZero);
% end
% P.fParfor = 0; %debug
try
    if P.fParfor %5x speed-up
        try
            parfor iSite = 1:nSites    
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                    spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);    
            end 
        catch
            P.fParfor=0;
            disp(lasterror);
        end
    end
    if ~P.fParfor
        for iSite = 1:nSites
            [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);            
        end
    end
catch
    disp(lasterr());
end

% merge parfor output and sort
viSpk = cell2mat_(cviSpkA);
vrSpk = cell2mat_(cvrSpkA);
viSite = cell2mat_(cviSiteA);
[viSpk, viSrt] = sort(viSpk);
vrSpk = vrSpk(viSrt);
viSite = viSite(viSrt);
end %func


function [viSpkA, vrSpkA, viSiteA] = spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P)
% spkLim = int32(abs(P.spkLim(1))*[-1,1]); %user shorter lim
spkLim = int32(abs(P.spkRefrac)*[-1,1]); %user shorter lim
% spkLim = int32(P.spkLim([1,end]));

vii1 = find(viSite == iSite);
viSpk1 = viSpk(vii1);
vrSpk1 = vrSpk(vii1);

% find neighbouring sites 
try     
    viSiteNear = P.miSites(:,iSite);
    vi2 = find(ismember(viSite, viSiteNear));
catch
    disp(lasterror());
    disp('spikeMerge: Using site numbers to determine neighbouring sites');
    siteLim = int32([-1,1]*P.maxSite); %int16 isn't faster   
    siteLim1 = iSite + siteLim;    
    vi2 = find(viSite>= siteLim1(1) & viSite <= siteLim1(2));
    %     vi2 = find(viSite == iSite | viSite == iSite-2 | viSite == iSite+2);
end
viSpk2 = viSpk(vi2); 
vrSpk2 = vrSpk(vi2);
viSite2 = viSite(vi2);
    
%     vlSpk2b(viSpk2/jitter) = 1;
vlSpk1 = false(size(viSpk1));
i2prev = 1;
for iSpk1=1:numel(viSpk1)
    iSpk11 = viSpk1(iSpk1);
    rSpk11 = vrSpk1(iSpk1); 

    % check for duplicate detection. search nearby spikes
    spkLim11 = iSpk11 + spkLim;
    [vii2, i2prev] = findRange_(viSpk2, spkLim11(1), spkLim11(2), i2prev); 
    [viSpk22, vrSpk22, viSite22] = select_vr(viSpk2, vrSpk2, viSite2, vii2);

    %invalid if larger (more negative) spike found
    if any(vrSpk22 < rSpk11), continue; end %wouldn't work for biopolar spike

    % check for equal amplitude, pick first occured
    vii22Eq = find(vrSpk22 == rSpk11);
    if numel(vii22Eq) > 1
        viSpk222 = viSpk22(vii22Eq);
        minTime = min(viSpk222);
        if minTime ~= iSpk11, continue; end 
        if sum(minTime == viSpk222) > 1 %pick lower site
            if min(viSite22(vii22Eq)) ~= iSite, continue; end
        end
    end
    vlSpk1(iSpk1) = 1; % set this spike as valid
end %for

% Trim
viiSpk1 = find(vlSpk1); %speed up since used multiple times
viSpkA = viSpk1(viiSpk1);
vrSpkA = vrSpk1(viiSpk1);
viSiteA = repmat(int32(iSite), size(viiSpk1));
% fprintf('.'); %progress. how about within parfor?
end %func


function [hPatch, hText] = plotProbe_(vcFile_prb)
if nargin<1, vcFile_prb='imec2.prb'; end
if matchFileExt(vcFile_prb, {'.bin', '.dat'})
    vcFile_prb = subsFileExt(vcFile_prb, '.prm');
end
if matchFileExt(vcFile_prb, '.prm')
    vcFile_prm = vcFile_prb;
    P = loadParam_(vcFile_prm);
    vcFile_prb = P.probe_file;
    if ~exist(vcFile_prb, 'file')
        vcFile_prb = replacePath(vcFile_prb, vcFile_prm);
    end
end

[viSite2Chan, mrSiteXY, vrSiteHW] = read_prb_file(vcFile_prb);

hFig = figure;
plot_probe(mrSiteXY, vrSiteHW, viSite2Chan);
vrPos0 = get(0, 'ScreenSize'); 
set(hFig, 'OuterPosition', vrPos0, 'Color', 'w');
axis equal;
set(hFig, 'Name', vcFile_prb, 'NumberTitle', 'off');
edit(vcFile_prb); %show probe file
figure(hFig);

end %func


function [Sevt, mrWav] = jrclust_event_(mrWav, varargin)
% all-in-one spike sorter. it assumes the binary file is already filtered
% and bad channels being excluded
% 
% vcFile: it can be filename, mrWav, or Sevt

% takes raw file
P = funcDefStr(funcInStr(varargin{:}), ...
    'spkLim', [-10 24], 'qqFactor', 5, 'qqSample', 8, ...
    'maxSite', 4, 'minAmpLimit', [], 'nChans', 120, ...
    'PERCENT', 0.2, 'nPcPerChan', 6, 'pca_n_waveforms_max', ...
    10000, 'dc', [], 'fCommonDc', 1, 'rho_target', 100, 'memLoad', 100, ...
    'fTwoSites', 0, 'vcCommonRef', 'median', 'vcFile', '');

if P.fPcaDetect
    fprintf('Extracting spike template ...'); t1 = tic;
    [P.vrPv, P.mrPv, P.vrPv_mean, mrWav_spk, Sevt1] = pca_detect_(mrWav, P);
    fprintf('\n\ttook %0.1fs.\n', toc(t1));
    P.spkThresh = median(Sevt1.vrThresh_site); % use global median for threshold
else
    P.mrPv = [];
    mrWav_spk = [];
end
Sevt = spikeDetect_Sevt_(mrWav, P);
Sevt.mrPv = P.mrPv;
if isfield(P, 'vrPv_mean'), Sevt.vrPv_mean = P.vrPv_mean; end
Sevt.mrWav_spk = mrWav_spk;
% Sevt = centroid_Sevt_(Sevt, mrWav, P); % spike centroid
% Sevt = getVpp_Sevt_(Sevt, mrWav, P); %needs viSite_spk (from centroid)
% Sevt = recenter_Sevt_(Sevt, P); %needs centroid
%  determine spike centroid

% Pricipal copmonent and apply to the waveform
% set2 serial run: 40.7s, parfor: 65.1s

% update spike location based on maximum site
% Sevt = Sevt_center_spike_(Sevt, mrWav, P);

fprintf('Getting feature...\n'); t3=tic;
Sevt = getFet_(Sevt, mrWav, P);
fprintf('Getting feature took %0.1fs\n', toc(t3));
Sevt.P = P; %save parameters
end %func


function Sevt = spikeDetect_Sevt_(mrWav, P)

P.spkThresh = P.spkThresh_uV / P.uV_per_bit;
[cviSpk, cvrSpk, vrThresh_site] = spikeDetect_(mrWav, P); 
[cviSpk, cvrSpk] = rejectSpk_Sevt_(cviSpk, cvrSpk, P); %reject spikes

% save spikes per site
if P.fSaveRawSpk
%     if strcmpi(P.vcSpatialFilter, 'none') && ~P.fPcaDetect
    Sevt.cviSpk_site = cviSpk;
    Sevt.cvrSpk_site = cvrSpk;
%     else
%         P1 = P;
%         P1.vrPv = []; %no matched filter
%         P1.vcSpatialFilter = 'none';
%         [Sevt.cviSpk_site, Sevt.cvrSpk_site, vrThresh_site] = ...
%             spikeDetect_(mrWav, P1); 
%     end
end
Sevt.vrThresh_site = vrThresh_site;
Sevt.vrThresh_uV = int16_to_uV_(vrThresh_site, P);

% combine spikes from multiple channels 
fprintf('Merging spikes...'); t2=tic;
[Sevt.viSpk, Sevt.vrSpk, Sevt.viSite] = spikeMerge_(cviSpk, cvrSpk, P);
fprintf('\t%d spikes found, took %0.1fs\n', numel(Sevt.viSpk), toc(t2));

% center the spike location at the amplitude max
% if ~strcmpi(P.vcSpatialFilter, 'none')
    Sevt = Sevt_center_spike_(Sevt, mrWav, P);
% end
end %func


function Sevt = Sevt_center_spike_(Sevt, mrWav, P)
% update max site and spike timing
% fUpdateTime = 0; fUpdateSite = ~strcmpi(P.vcSpatialFilter, 'none');
return; %dont do it
fUpdateSite = 1; fUpdateTime = 1;

fprintf('Centering spikes\n'); t1=tic;
viSite0 = Sevt.viSite; %before
it0 = 1 - P.spkLim_fet(1);
for iSite=1:size(mrWav,2)
    viiSpk1 = find(viSite0==iSite);
    if isempty(viiSpk1), continue; end
    viTime1 = Sevt.viSpk(viiSpk1); 
	viSite1 = P.miSites(:,iSite);
    trWav1 = mr2tr3(mrWav, P.spkLim_fet, viTime1, viSite1);
%     trWav1 = cumsum(trWav1,1)/4;
    mrVpp1 = shiftdim(max(trWav1,[],1)-min(trWav1,[],1));
%     mrVpp1 = shiftdim(abs(min(trWav1,[],1)));
%     mrVpp1 = abs(shiftdim(max(trWav1,[],1) - min(trWav1,[],1)));
%     mrVpp1 = abs(squeeze(max(trWav1(it0+1:end,:,:)) - min(trWav1(1:it0,:,:))));
%     [vrSpk1, viiSite1] = min(mrWav(viTime1, viSite1), [], 2);
%     if iscolumn(mrVpp1), mrVpp1 = mrVpp1'; end
    
    % update max site
    [Sevt.vrSpk(viiSpk1), viiSite1] = max(mrVpp1, [], 2);
%     Sevt.vrSpk(viiSpk1) = abs(vrSpk1); %spike amplitude
    if fUpdateSite
        Sevt.viSite(viiSpk1) = viSite1(viiSite1);
    end
    
    % update timing
    if fUpdateTime
        dimm1 = size(trWav1);
        trWav1 = reshape(trWav1, dimm1(1), []);
        mrWav1 = trWav1(:, sub2ind(dimm1(2:3), 1:dimm1(2), viiSite1(:)'));
%         mrWav1 = cumsum(mrWav1,1); % align at the integrated minimum
        [~, viMin1] = min(mrWav1);
%         [~, viMax1] = max(rms_filt(single(mrWav1), P.nRms_filt));
        Sevt.viSpk(viiSpk1) = viTime1 + int32(viMin1' + P.spkLim_fet(1)-1);
    end
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [cviSpk, cvrSpk] = rejectSpk_Sevt_(cviSpk, cvrSpk, P)
% reject spikes based on the reference signal
if isempty(P.rejectSpk_mean_thresh), return; end
if P.rejectSpk_mean_thresh <= 0, return; end
try
    vlKeep = load_ref_(subsFileExt(P.vcFile_prm, '.ref'), P);

%     if ~isempty(vlKeep)
    for iSite=1:numel(cviSpk)
        viSpk1 = cviSpk{iSite};
        if isempty(viSpk1), continue; end %skip if no spikes found
        [cviSpk{iSite}, cvrSpk{iSite}] = ...
            select_vr(viSpk1, cvrSpk{iSite}, find(vlKeep(viSpk1)));
    end
%     end
catch
    disp(['load_ref: ', lasterr()]);
end
end % func


function [csFile_merge, vcDir, csFile_merge1] = dir_file_(vcFile_dir, fSortByDate)
% search for files and sort by date
if nargin<2, fSortByDate=1; end

[vcDir, ~, ~] = fileparts(vcFile_dir);
vsDir = dir(vcFile_dir);
vrDatenum = cell2mat({vsDir.datenum});
csFile_merge = {vsDir.name};
if fSortByDate
    [~,ix] = sort(vrDatenum, 'ascend');
    csFile_merge = csFile_merge(ix);
end
csFile_merge1 = csFile_merge; 
csFile_merge = cellfun(@(vc)[vcDir, filesep(), vc], csFile_merge, 'UniformOutput', 0);
end %func


function merge_binfile_(vcFile, csFile_merge)
% skip if file already exists
if exist(vcFile, 'file'), 
    disp([vcFile ' already exists']); 
    return; 
end

% get a list of flies
if ~iscell(csFile_merge)
    fSortByDate = 1;
    [csFile_merge1, vcDir, csFile_merge] = dir_file_(csFile_merge, fSortByDate);
    fSuccess = merge_binfile_system_(vcFile, csFile_merge, vcDir);
else    
    csFile_merge1 = csFile_merge;
    fSuccess = merge_binfile_system_(vcFile, csFile_merge);
end

% merge files
if ~fSuccess %problem with merging
    merge_binfile_slow_(vcFile, csFile_merge1);
end

% copy meta file
try
    vcFile1_meta = subsFileExt(csFile_merge1{1}, '.meta');
    vcFile_meta = subsFileExt(vcFile, '.meta');
    copyfile(vcFile1_meta, vcFile_meta, 'f');
catch
    disp('meta file does not exist for merging');
end

% Merge LFP file for IMEC3 probe
try
    if ~isempty(strfind(lower(vcFile), '.imec.ap.bin'))    
        func_ap2lf = @(x)strrep(lower(x), '.imec.ap.bin', '.imec.lf.bin');
        vcFile_lf = func_ap2lf(vcFile);
        csFile_merge_lf = cellfun(@(x)func_ap2lf(x), csFile_merge1, 'UniformOutput', 0);
        merge_binfile_(vcFile_lf, csFile_merge_lf);
    end
catch
    disp('Merge LFP file error for IMEC3.');
end
end %func


function fSuccess = merge_binfile_system_(vcFile, csFile_merge, vcDir)
% Merge files (csFile_merge) to vcFile
if nargin<3, vcDir = ''; end
if ~isempty(vcDir)
    vcDir_prev = pwd();
    cd(vcDir);
else
    vcDir_prev = '';
end

vcCmd = 'copy /b ';
for iFile = 1:numel(csFile_merge)
    if iFile < numel(csFile_merge)
        vcCmd = sprintf('%s"%s"+', vcCmd, csFile_merge{iFile});
    else
        vcCmd = sprintf('%s"%s" "%s"', vcCmd, csFile_merge{iFile}, vcFile);
    end
end
try
    t1 = tic;
    status = system(vcCmd);
    fprintf('\tFile merging took %0.1fs\n', toc(t1));
    cd(vcDir_prev);
    fSuccess = status == 0;
catch
    fSuccess = 0;
end
if ~isempty(vcDir_prev), cd(vcDir_prev); end
end %func


function merge_binfile_slow_(vcFile, csFile_merge)
% system call is 3.1x faster but this has no limits on number of files
% merged

t1 = tic;
fidw = fopen(vcFile, 'W');    
for iFile=1:numel(csFile_merge)
    vcFile1 = csFile_merge{iFile};
    fidr = fopen(vcFile1, 'r');
    vr = fread(fidr, inf, 'int16=>int16');
    fclose(fidr);    
    fwrite(fidw, vr, 'int16');
    disp(vcFile1);
end
fclose(fidw);
fprintf('\tFile merging took %0.1fs\n', toc(t1));
end %func


function [mrWav, mrLfp, mrAux, viSiteZero] = export_spk_(vcFin, varargin)
% export to .spk format
% global mrWav
P = funcDefStr(funcInStr(varargin{:}), 'sRateHz', 25000, ...
    'nChans', 129, 'freqLim', [300 3000], 'filtOrder', 3, 'fEllip', 1, ...
    'freqLimNotch', [], 'freqLimStop', [], ...
    'vcCommonRef', 'median', 'fCheckSites', 1, 'maxLfpSdZ', 4.5, ...
    'nLfpSubsample', 10, 'fSaveSpk', 1, 'probe_type', 'imecii', ...
    'viSite2Chan', [], 'mrSiteXY', [], 'vrSiteWH', [], 'viSiteZero', [], ...
    'uV_per_bit', (4/2^16/200*1e6), 'gain_boost', 1);

t1=tic;
[vcDir, vcFile, vcExt] = fileparts(vcFin);
if isempty(vcDir), vcDir = '.'; end

if ~isempty(P.csFile_merge)
    % create a one big .bin file
    merge_binfile_(vcFin, P.csFile_merge);
end
if ~exist(vcFin, 'file'), error(['File not found: ' vcFin]); end

% filter and load 
[mrWav, mrLfp, mrAux, P, vrWav_ref] = load_bin_(vcFin, P);

% reject global outlier
% viRefOut = vrWav_ref>P.
% for iChan=1:size(mrWav,2)
%     mrWav(iChan,viRefOut) = 0;
% end


% Set bad sites to zero
% if ~isempty(P.viSiteZero)
viSiteZero = P.viSiteZero;
% elseif P.fCheckSites
%     viSiteZero = detectBadSites_(mrLfp, P);
% else
%     viSiteZero = [];
% end
% if ~isempty(viSiteZero)
%     mrWav(:, viSiteZero) = 0;
%     if ~isempty(mrLfp), mrLfp(:, viSiteZero) = 0; end
%     disp(['Sites ', sprintf('%d, ', viSiteZero), 'set to zero.']);
% end
% P.viSiteZero = viSiteZero;

% determine reference channel
% nChans = size(mrWav,2);
% viChan = setdiff(1:nChans, P.viSiteZero);
% if P.fCheckSites && ~isempty(viSiteZero)
%     fprintf('Averaging across channel '); t2=tic;    
%     vrWav_ref = [];
%     for iCh = viChan
%         if isempty(vrWav_ref)
%             vrWav_ref = single(mrWav(:,iCh)); %global mean
%         else
%             vrWav_ref = vrWav_ref + single(mrWav(:,iCh)); %global mean
%         end
%         fprintf('.');
%     end
%     vrWav_ref = vrWav_ref / numel(viChan);
%     fprintf(' took %0.1fs\n', toc(t2));    
% end

% filter LFP
% if ~isempty(P.freqLim_lfp) || ~isempty(P.freqLimNotch_lfp)
%     fprintf('Filtering LFP '); t2=tic;
%     for iCh = viChan
%         if isempty(mrLfp), break; end
%         vrLfp1 = single(mrLfp(:,iCh));
%         vrLfp1 = firfilt(vrLfp1, P.freqLim_lfp, P.sRateHz_lfp);    
%         vrLfp1 = notchfilt(vrLfp1, P.freqLimNotch_lfp, P.sRateHz_lfp);
%         mrLfp(:,iCh) = int16(vrLfp1);
%         fprintf('.');
%     end
%     fprintf(' took %0.1fs\n', toc(t2)); 
% end

if strcmpi(P.vcDataType, 'single')
    mrWav = int16(mrWav / P.uV_per_bit);
    vrWav_ref = int16(vrWav_ref / P.uV_per_bit);
else
    vrWav_ref = int16(vrWav_ref);
end
if P.fSaveSpk      
    write_bin_(subsFileExt(P.vcFile_prm, '.spk'), mrWav);
    write_bin_(subsFileExt(P.vcFile_prm, '.lfp'), mrLfp);
    write_bin_(subsFileExt(P.vcFile_prm, '.ref'), vrWav_ref);
    write_bin_(subsFileExt(P.vcFile_prm, '.aux'), mrAux);
end
fprintf('Exported %s (took %0.1fs)\n', vcFin, toc(t1));
end %func


function write_bin_(vcFile, mr)
% non-blocking write
% start(timer('TimerFcn', @(e,h)write_bin_timer(vcFile, mr), 'StartDelay', .1));
write_bin_timer(vcFile, mr);
end %func


function viSiteZero = detectBadSites_(mrLfp, P)
if isempty(mrLfp), viSiteZero = []; return; end

% mrLfp = mrWav(1:P.nLfpSubsample:end,:);
vrLfp_mean = int16(mean(mrLfp,2));
vrLfpSdZ = zscore(std(bsxfun(@minus, mrLfp, vrLfp_mean))); %remove mean across chan    
viSiteZero = find(vrLfpSdZ > P.maxLfpSdZ);
end %func


function benchtest_(vcFile_prm)
% plots RMS of the raw band
% workspace export RMS
% Bench test result

% if matchFileExt(vcFile_prm, {'.bin', '.dat'})
%     vcFile_prm = subsFileExt(vcFile_prm, '.prm');
% end
    
P = loadParam_(vcFile_prm);

mrWav = loadWavFile(P.vcFile, P.nChans, inf, P.vcDataType, P.viSite2Chan);

vrRms = zeros(size(mrWav,2), 1, 'single');
for iCh=1:size(mrWav,2)
    vrRms(iCh) = std(single(mrWav(:,iCh))) * P.uV_per_bit;
end
bitDepth = mr_bitDepth(mrWav);
if bitDepth > 10
    %Bill's spikeglx software number representation changed
    vrRms = vrRms / 2^6; 
end

uVrms = vrRms;
printstat(uVrms);

uVpp = vrRms * 2 * sqrt(2);
printstat(uVpp);

figure; stem(vrRms); xlabel('Site #'); ylabel('Vrms (uV)');
set(gcf, 'Name', P.vcFile);

assignWorkspace(vrRms, mrWav);

end %func


% function vcFile_prb = guess_prb_file_(vcFile_bin)
% if matchFileExt(vcFile_bin, '.imec.bin')
%     
% else
%     
% end
% end %func


function sitestat_(vcFile_prm)
vcFile_evt = subsFileExt(vcFile_prm, '_evt.mat');
Sevt = load_evt_(vcFile_prm);
viValid_site = find(cellfun(@(x)~isempty(x), Sevt.cviSpk_site));
try
    tDur = Sevt.P.S_imec3.fileTimeSecs; %for imec3 only
catch
    tDur = max(cellfun(@(x)max(x), Sevt.cviSpk_site(viValid_site))) / Sevt.P.sRateHz;
end
% counts_site = cellfun(@(x)numel(x), Sevt.cvrSpk_site(viValid_site));
rate = counts_site / tDur;
printstat(rate);
viSites_plot = viValid_site; %Sevt.P.viSite2Chan(viValid_site);
amp90 = cellfun(@(x)quantile(abs(x), .9), Sevt.cvrSpk_site(viValid_site));
amp90 = int16_to_uV_(amp90, Sevt.P);
printstat(amp90);

try
    rms_qq = Sevt.vrThresh_uV(viValid_site) / Sevt.P.qqFactor;
    SNR = amp90 ./ rms_qq;
    printstat(SNR);
    printstat(rms_qq);
catch
    ;
end

figure; 
set(gcf, 'Name', vcFile_prm, 'Color', 'w');
ax=[];
ax(1) = subplot(211); 
bar(viSites_plot, amp90,1); ylabel('90th ampl (uV)'); xlabel('Sites'); grid on; set(gca,'YScale','log');
ax(2) = subplot(212); 
bar(viSites_plot, counts_site, 1); ylabel('Spike counts'); xlabel('Sites'); grid on; set(gca,'YScale','log');
linkaxes(ax,'x');
xlim([min(viSites_plot), max(viSites_plot)]);
end %func


function [vii, ic] = findRange_(vi, a, b, i)
% a: start value
% b: end value
% i: index of vi to start searching
% vii: index of vi that is between a and b
% vii = [];
ic = 1;
% try
n = numel(vi);
ia = 0;
while 1
    if vi(i) >= a && ia == 0, ia = i; end      
    if vi(i) < a, ic = i; end    
    i = i + 1;    
    if i > n, ib = n; break; end    
    if vi(i) > b, ib = i-1; break; end    
end
if ia > 0
    vii = ia:ib;
else
    vii = [];
end
end %func


function viData = readChan_whisper_(vcFile, nChans, iChan)
% iChan starts with 0
Sfile = dir(vcFile);
nSamples = floor(Sfile(1).bytes / nChans / 2);
fid = fopen(vcFile);
fseek(fid, iChan * 2, 'bof'); %go to position
viData = fread(fid, nSamples, '*int16', 2 * nChans-2);
fclose(fid);
end %func


function nBytes = file_bytes_(vcFile)
Sfile = dir(vcFile);
if isempty(Sfile)
    nBytes = 0;
else
    nBytes = Sfile(1).bytes;
end
end


function viMid = find_pulse_mid_(viSync)
thresh_sync = mean([max(viSync), min(viSync)]);
viSync1 = diff(viSync > thresh_sync);
viUp = find(viSync1 > 0);
viDn = find(viSync1 < 0);
if viDn(1) < viUp(1), viDn(1) = []; end
if viDn(end) < viUp(end), viUp(end) = []; end
viMid = round((viUp + viDn) / 2 + 1);
end %func


function led_sync_(P)
if ~exist(P.vcFile_vid, 'file')
    error(['video file does not exist: ', P.vcFile_vid]); 
end

vidobj = VideoReader(P.vcFile_vid);

% get a crop
trImg1 = readFrame(vidobj);
hFig = create_figure('iFig', 3001, 'name', ['Probe map; ', P.vcFile_vid], 'pos', [0 0 .5 1], 'fMenubar', 1, 'fToolbar', 1);
imshow(trImg1);
uiwait(msgbox('draw a rectangle around LED'));
hPoly = imrect;
mrPolyPos = getPosition(hPoly);
xlim = round(cumsum(mrPolyPos([1,3]))); 
ylim = round(cumsum(mrPolyPos([2,4]))); 
viY = ylim(1):ylim(2);
viX = xlim(1):xlim(2);
mr1 = trImg1(viY, viX, 1);

iFrame  = 1;
nFrames = ceil(vidobj.Duration * vidobj.FrameRate);
tnLed = zeros([size(mr1), nFrames], 'like', mr1);
tnLed(:,:,iFrame) = mr1;
% vrT = (1:nFrames) / vidobj.FrameRate;
hTitle = title(sprintf('Frame 1, 0s, 0.0%%'));

% collect led
while hasFrame(vidobj)
    tnImg1 = readFrame(vidobj);
    tnLed(:,:,iFrame) = sum(tnImg1(viY, viX), 3);    
    iFrame = iFrame + 1;
    if mod(iFrame, vidobj.FrameRate)==0
        set(hTitle, 'String', sprintf('Frame %d, %0.1fs, %0.1f%%', iFrame, iFrame / vidobj.FrameRate, iFrame/nFrames*100));
        drawnow; 
    end
end
tnLed = tnLed(:,:,1:iFrame-1); %trim unused

viSync = readChan_whisper_(P.vcFile, P.nChans, P.iChan_vid);
vrTime_sync = find_pulse_mid_(viSync) / P.sRateHz;
assignWorkspace(tnLed, viSync, vrTime_sync);
implay(tnLed);

% viInt_ephys = round(diff(vrTime_sync));
% viFrame = [3 84 144 214 284 364 434 495 574 634 704 834 884 955 1024 1105 1165 1225 1294];
% fprintf('%d', round(diff(viFrame)/10)); fprintf('\n');
% findstr(sprintf('%d', viInt_ephys), '5778667')
% 12th video sequence lines with 3rd interval

end %func


function showTraj_(P)
% plot trajectory and use bonsai
vidobj = VideoReader(P.vcFile_vid);
img = readFrame(vidobj);
Sxy = load(P.vcFile_bonsai);
vrX = Sxy.animalpos(:,2);
vrY = Sxy.animalpos(:,3);

figure; imshow(img);
hold on; plot(vrX, vrY);
title(P.vcFile_vid);
set(gcf, 'Name', P.vcFile_prm);

assignWorkspace(vrX, vrY, P);
end %func


function trFet = fet_diff248_(mrWav, viSpk, P)
vrWav0 = mrWav(viSpk, :);
nT = size(mrWav,1);

viSpk8 = viSpk + 8;
viSpk4 = viSpk + 4;
viSpk1 = viSpk + 2;
viSpk8(viSpk8>nT) = nT;
viSpk4(viSpk4>nT) = nT;
viSpk1(viSpk1>nT) = nT;

trFet = cat(3, mrWav(viSpk1,:)-vrWav0, mrWav(viSpk4,:)-vrWav0, mrWav(viSpk8,:)-vrWav0);
trFet = single(permute(trFet, [3,2,1]));
end %func


function commit_version_()
% update jrclust_alpha to jrclust directory
t1 = tic;
S_cfg = read_cfg_();
csCopy = S_cfg.sync_list;
vcDest = S_cfg.path_dropbox;

if ~strcmpi(pwd(), S_cfg.path_alpha), disp('must commit from alpha'); return; end

disp(['Commiting to ', vcDest]);

% delete empty files and temp files
delete_files_(find_empty_files_());
% delete('temp_jrclust_*.m');

for iCopy = 1:numel(csCopy);
    vcEval1 = sprintf('copyfile ''%s'' ''%s'' f;', csCopy{iCopy}, vcDest);    
    try        
        eval(vcEval1);
        fprintf('\t%s\n', vcEval1);  
    catch
        fprintf(2, '\tError copying %s\n', csCopy{iCopy});
    end
end
edit change_log.txt
fprintf('Commited, took %0.1fs.\n', toc(t1));
end %func


function csFiles = find_empty_files_(vcDir)
% find files with 0 bytes
if nargin==0, vcDir = pwd(); end
vS_dir = dir(vcDir);
viFile = find([vS_dir.bytes] == 0 & ~[vS_dir.isdir]);
csFiles = {vS_dir(viFile).name};
csFiles = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


function delete_files_(csFiles)
% if ischar(csFiles)
%     vS_dir = dir(csFiles);
%     vS_dir = vS_dir(~(vS_dir.isdir));
%     csFiles = cellfun(@(vc)[vcDir, filesep(), vc], {vS_dir.name}, 'UniformOutput', 0);
% end
for iFile = 1:numel(csFiles)
    try
        delete(csFiles{iFile});
    catch
        ;
    end
end
end %func


function update_version_(S_cfg)
% update jrclust_alpha to jrclust directory
% vcSource = 'C:\Dropbox (HHMI)\Git\jrclust';
if nargin<2, S_cfg = read_cfg_(); end
vcSource = S_cfg.path_dropbox;
if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
vcBackup = S_cfg.path_backup;
mkdir(vcBackup);

t1 = tic;
fprintf('Updating from %s...\n', vcSource);
csCopy = S_cfg.sync_list; 
for iCopy = 1:numel(csCopy);
    vcEval1 = sprintf('copyfile ''%s\\%s'' .\\ f;', vcSource, csCopy{iCopy});
    try
        try eval(sprintf('copyfile ''.\\%s'' ''%s\\'' f;', csCopy{iCopy}, vcBackup)); catch; end
        eval(vcEval1); 
        fprintf('\t%s\n', vcEval1);  
    catch
        fprintf(2, '\tError updating %s\n', csCopy{iCopy});  
    end
end
try
    compile_cuda_(S_cfg);
catch
    fprintf(2, 'CUDA code compilation error.\n');
end
fprintf('Updated, took %0.1fs.', toc(t1));
fprintf('\tPrevious files backed up to %s\n', vcBackup);
edit change_log.txt
end %func


function plot_event_rate_(vcFile_prm)
% event-rate command
P = file2struct(vcFile_prm);
Sevt = load_evt_(P.vcFile_prm);

viSite=double(Sevt.viSite); 
viSite_site = 1:max(viSite);
vnCnt_site = hist(viSite, viSite_site); 
t_dur = double(max(Sevt.viSpk)) / P.sRateHz;
vrRate_site = vnCnt_site / t_dur;
figure;  bar(viSite_site, vrRate_site, 1); xlabel('Site #'); ylabel('Event Rate (Hz)'); set(gca, 'YScale','log'); grid on;
set(gcf, 'Name', P.vcFile, 'Color', 'w');
S_site = makeStruct(viSite_site, vrRate_site);
assignWorkspace(S_site, Sevt);
end %func


function plot_activity_(P)
% plot activity as a function of depth and time
tbin = 10; %activity every 10 sec
% plot activity as a function of time
% vcFile_evt = subsFileExt(P.vcFile_prm, '_evt.mat');
Sevt = load_evt_(P.vcFile_prm);
nSites = numel(Sevt.cvrSpk_site);
tdur = max(cell2mat(cellfun(@(x)double(max(x)), Sevt.cviSpk_site, 'UniformOutput', 0))) / P.sRateHz;
nTime = ceil(tdur / tbin);


mrAmp90 = zeros(nTime, nSites);
lim0 = [1, tbin * P.sRateHz];
for iSite=1:nSites
    vrSpk1 = Sevt.cvrSpk_site{iSite};    
    if isempty(vrSpk1), continue; end
    viSpk1 = Sevt.cviSpk_site{iSite};
    for iTime=1:nTime
        lim1 = lim0 + (iTime-1) * lim0(2);
        vrSpk11 = vrSpk1(viSpk1 >= lim1(1) & viSpk1 <= lim1(2));
        if isempty(vrSpk11),  continue; end
        mrAmp90(iTime, iSite) = quantile(abs(vrSpk11), .9);
    end
end %for
mrAmp90=mrAmp90';

vlSite_left = P.mrSiteXY(:,1) == 0;
vrSiteY = P.mrSiteXY(:,2);
figure; 
subplot 121; imagesc(mrAmp90(vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites');
subplot 122; imagesc(mrAmp90(~vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(~vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites');

[~, iSite_center] = max(mean(mrAmp90,2));
viSiteA = iSite_center + [-2:2]; %look neighbors
mrAmp90a = mrAmp90(viSiteA, :);
vrCentroid = bsxfun(@rdivide, sum(bsxfun(@times, mrAmp90a.^2, vrSiteY(viSiteA))), sum(mrAmp90a.^2));
hold on; plot((1:nTime) * tbin, vrCentroid);

end %func


function plot_activity1_(P)
tbin = 5; %activity every 10 sec
% plot activity as a function of time
% vcFile_evt = subsFileExt(P.vcFile_prm, '_evt.mat');
Sevt = load_evt_(P.vcFile_prm);
nSites = numel(Sevt.cvrSpk_site);
tdur = max(cell2mat(cellfun(@(x)double(max(x)), Sevt.cviSpk_site, 'UniformOutput', 0))) / P.sRateHz;
nTime = ceil(tdur / tbin);

% mrAmp90 = zeros(nTime, nSites);
lim0 = [1, tbin * P.sRateHz];
vrCentroid = zeros(nTime, 1);
vrSiteY = P.mrSiteXY(:,2);
for iTime=1:nTime
    lim1 = lim0 + (iTime-1) * lim0(2);
    vlSpk1 = Sevt.viSpk >= lim1(1) & Sevt.viSpk <= lim1(2);
    vrSpk1 = single(abs(Sevt.vrSpk(vlSpk1)));
    viSite1 = Sevt.viSite(vlSpk1);
    vlSpk11 = vrSpk1 > quantile(vrSpk1, .9);
    vrSpk11 = vrSpk1(vlSpk11);
    viSite11 = viSite1(vlSpk11);
    vrCentroid(iTime) = sum(vrSpk11 .* vrSiteY(viSite11)) ./ sum(vrSpk11);
end

figure; plot((1:nTime) * tbin, vrCentroid);

end %func


function [mnWav, mnLfp, mnAux, P, vrRef] = load_bin_(vcFile_bin, P)
% filter is done while loading
LOAD_FACTOR = 8; %GPU memory usage factor. 4x means 1/4 of GPU memory can be loaded
fImec3 = match_postfix(vcFile_bin, '.imec.ap.bin') || match_postfix(vcFile_bin, '.imec.bin') || P.nChans >= 276;
if fImec3, P.nSkip_lfp = 12; end %extra protection. machine set parameter
 
P.nSamples_load = floor(P.MAX_BYTES_LOAD / bytePerData(P.vcDataType) / P.nChans / P.nSkip_lfp) * P.nSkip_lfp;
if isempty(P.nSamples_load), P.nSamples_load = floor(mem_max_(P) / LOAD_FACTOR); end
P.nSamples_load = 2^(nextpow2(P.nSamples_load)-1); %GPU friendly for FFT operation

if fImec3
    %imec 3 SPIKEGLX format
%     P.nSamples_load = 30000 * 1 * 12;
    [mnWav, mnLfp, mnAux, Smeta, vrRef] = import_imec3(vcFile_bin, P); %load_imec3 is more memory intensive
    P.S_imec3 = Smeta.S_imec3; % Smeta is loaded again as S_imec3. parse imroTbl   
else
    % imec phase II and other whisper format
    P.viChan = P.viSite2Chan;  
    P.viChan_lfp = P.viChan;
%     P.nSamples_load = 25000 * 2 * 10;
    [mnWav, mnLfp, mnAux, vrRef] = import_whisper(vcFile_bin, P);
end

% if strcmpi(P.vcDataType, 'uint16') 
%     % convert data type to int16
%     mnWav = mnWav - 2^15; 
%     mnLfp = mnLfp - 2^15;
%     vrRef = vrRef - 2^15;
% end
end


function mnLfp = subsample_lfp_(mnWav, nSkip)   
nSamples = size(mnWav,1);
nSamples_lfp = floor(nSamples / nSkip);
viLfp = (0:nSamples_lfp-1) * nSkip + 1;
mnLfp = mnWav(viLfp,:);
end %func


function save_struct_(vcFile, S)
% Write a stuct S to a file
% none-blocking
% fprintf('Saving to %s...', vcFile); 
% t1=tic;

% start(timer('write_bin_timer
% save(vcFile, '-struct', 'S', '-v7.3');

start(timer('TimerFcn', @(e,h)write_bin_timer(vcFile, S), 'StartDelay', .1));


% fprintf(' took %0.1fs\n', toc(t1));
end %func


function vr_uV = int16_to_uV_(vn_int16, P)
try
    % imec3 format
    vr_uV = single(vn_int16) .* P.S_imec3.vrScale_ap;
catch
    vr_uV = single(vn_int16) * P.uV_per_bit;
end
end %func


function site_offset_(vcFile_prm)
vcFile_bin = file_prm2bin_(vcFile_prm);

% read the 1-3 min data
[mnWav, mnLfp, mnAux, Smeta] = import_imec3(vcFile_bin, 'tlim_load', [60, 180]);
assignWorkspace(mnWav, mnLfp, mnAux, Smeta);
end %func


function importTraces_ap_(vcFile_prm)
vcFile_bin = file_prm2bin_(vcFile_prm);

% read the 1-3 min data
[mnWav, mnLfp, mnAux, Smeta] = load_bin_(vcFile_bin);
assignWorkspace(mnWav, mnLfp, mnAux, Smeta);
end %func



function load_raw_(vcFile_prm)
vcFile_bin = file_prm2bin_(vcFile_prm);

% read the 1-3 min data

[mnWav, mnLfp, mnAux, Smeta] = import_imec3(vcFile_bin);

assignWorkspace(mnWav, mnLfp, mnAux, Smeta);

end %func


function plotLfpCorr_(vcFile_prm)
vcFile_bin = file_prm2bin_(vcFile_prm);

% read the 1-3 min data
P = loadParam_(vcFile_prm);
[mnWav, mnLfp, mnAux, Smeta] = import_imec3(vcFile_bin, 'tlim_load', [60, 180]);

mrLfp = zeros(size(mnLfp,1), max(Smeta.viSites), 'single');
mrLfp(:, Smeta.viSites) = mnLfp;
% P.sRateHz_lfp = 30000/12;
mrCoh = mt_cpsd_mr(mrLfp, 'sRateHz', 30000/12, 'freqLim', P.freqLim_corr, 'freqLim_excl', [58 62]);

% mrLfp = bandfilt(single(mrLfp), 'freqLim', P.freqLim_corr, 'sRateHz', Smeta.sRateHz/P.nSkip_lfp);
% figure; imagesc(corr(single(mnLfp)));
figure; imagesc(mrCoh);
title(sprintf('LFP coherence (%0.1f-%0.1f Hz)', P.freqLim_corr));
set(gcf, 'Name', vcFile_bin, 'Color', 'w');
axis xy;
xlabel('Site #');
ylabel('Site #');
colormap jet; caxis([0 1]);
axis square;

end %func


function [vcFile_bin, P] = file_prm2bin_(vcFile_prm)
if matchFileExt(vcFile_prm, '.prm')
    P = loadParam_(vcFile_prm);
    vcFile_bin = P.vcFile;
else
    vcFile_prm = subsFileExt(vcFile_prm, '.prm');
    [vcFile_bin, P] = file_prm2bin_(vcFile_prm);    
end
end %func


function [vcFile_lfp, P] = imec_sites_(P)
vcFile_lfp = lower(P.vcFile);
viSites_dat = [];
% mrSingleColumn = [];
% return site2chan map and reference sites
if ~isempty(strfind(vcFile_lfp, '.imec.bin'))
    % phase 3
    nSites = floor(P.nChans/2);
    viSite2Chan = (1:nSites) + nSites;
    viChan_ref = [37 76 113 152 189 228 265 304 341 380];
    viChan_ref(viChan_ref>nSites) = [];
    nSkip_lfp = 12;
    nSkip_site = 4;
elseif ~isempty(strfind(vcFile_lfp, '.imec.ap.bin'))
    nSites = P.nChans - 1;
    viSite2Chan = 1:nSites;
    viChan_ref = [37 76 113 152 189 228 265 304 341 380];
    viChan_ref(viChan_ref>nSites) = [];
    vcFile_lfp = strrep(vcFile_lfp, '.imec.ap.bin', '.imec.lf.bin');
    nSkip_lfp = 1;
    nSkip_site = 4;
elseif ~isempty(strfind(vcFile_lfp, '.imec_lfp.bin'))
    % Nick's LFP format
    nSites = P.nChans;
    viSite2Chan = 1:nSites;
    viChan_ref = [37 76 113 152 189 228 265 304 341 380];
    viChan_ref(viChan_ref>nSites) = [];
    nSkip_lfp = 1;
    nSkip_site = 2; %2 for even/odd and 4 for four-column
elseif ~isempty(strfind(vcFile_lfp, '.imec_lf.bin'))
    % My LFP format. 2500 sampling rate
    nSites = P.nChans - 1;
    viSite2Chan = 1:nSites;
    viChan_ref = [37 76 113 152 189 228 265 304 341 380];
    viChan_ref(viChan_ref>nSites) = [];
    nSkip_lfp = 1;
    nSkip_site = 2;    
else
    %assume phase 2
    % phase 2
    viSite2Chan = 1 + [ ...
        40  103 39	104	41	102	38	105	42	101	37	106	43	100	36	107 ...
        44	99	35	108	45	98	34	109	46	97	33	110	47	96	32	111 ...
        48	127	63	112	49	126	62	113	50	125	61	114	51	124	60	115	...
        52	116	59	123	53	117	58	122	54	118	57	121	55	119	56	120	...
        8	71	7	72	9	70	6	73	10	74	5	69	11	75	4	68 ...	
        12	76	3	67	13	77	2	66	14	78	1	65	15	79	0	64 ...	
        16	80	31	95	17	81	30	94	18	82	29	93	19	83	28	92 ...	
        20	84	27	91	21	85	26	90	22	86	25	89	23	87	24	88];
    viChan_ref = [1 18 33 50 65 82 97 114];    
    nSkip_lfp = P.nSkip_lfp;
    if matchFileExt(P.vcFile, '.dat') % Joao's format
        viSites_dat = setdiff(1:numel(P.viSite2Chan), P.viChan_ref);
    end    
    nSkip_site = 2;
end

P.viSite2Chan = viSite2Chan;
P.viChan_ref = viChan_ref;
P.nSkip_lfp = nSkip_lfp;
P.viSites_dat = viSites_dat;
P.viSite_odd = 1:2:numel(P.viSite2Chan);
P.viSite_even = 2:2:numel(P.viSite2Chan);
P.cviSite_track = partition_sites_(numel(P.viSite2Chan), nSkip_site);
end %func


function cviSite = partition_sites_(nSites, nSkip)
cviSite = cell(nSkip, 1);
for iSkip=1:nSkip
    cviSite{iSkip} = iSkip:nSkip:nSites;
end
end %func


function [mr1, viChan_ref_minus, viChan_ref_plus] = interp_ref_(mr1, viChan_ref)
nSites = size(mr1,2);
[viChan_ref_minus, viChan_ref_plus] = neighbor_sites_(nSites, viChan_ref);
mr1(:,viChan_ref) = (mr1(:,viChan_ref_minus) + mr1(:,viChan_ref_plus))/2;
end %func


function [viChan_ref_minus, viChan_ref_plus, viChan_ref] = neighbor_sites_(nSites, viChan_ref)
if nargin<2, viChan_ref = 1:nSites; end
viChan_ref_minus = viChan_ref - 2;
viChan_ref_plus = viChan_ref + 2;
viChan_ref_minus(viChan_ref_minus<1) = viChan_ref(viChan_ref_minus<1) + 2;
viChan_ref_plus(viChan_ref_plus>nSites) = viChan_ref(viChan_ref_plus>nSites) - 2;
end %func


function mr_column = compute_column_(mr, cviSite_track)
nSites_col = numel(cviSite_track{1});
if nSites_col == 64 
    %phase 2
%     mrSingleColumn = ones(nSites_col, 2); 
%     mr_column = mr(:, cviSite_track{1}) + mr(:, cviSite_track{2});
    mr_column = mr(:, cviSite_track{1}); %use the best edge
%     mr_column = 1./ (1./mr(:, cviSite_track{1}) + 1./mr(:, cviSite_track{2})); %v~1/r^2
%     mr_column = (mr(:, cviSite_track{1}).^-2 + mr(:, cviSite_track{2}).^-2) .^ 2; %v~1/r    
%     V1 = double((mr(:, cviSite_track{1})));
%     V2 = double((mr(:, cviSite_track{2})));
%     vl_zero_1 = V1==0;
%     vl_zero_2 = V2==0;
%     hmean = 2./(1./V1+1./V2);
    
    
%     amean = (V1+V2)/2;
%     min1 = min(min(V1(:)), min(V2(:))) + 1;
%     V1 = V1 + min1;
%     V2 = V2 + min1;
%     hmean = 2./(1./V1+1./V2) - min1;
%     mr_column = hmean;
% 
%     figure; plot(hmean(:,10)); hold on; plot(amean(:,10)); plot(V1(:,10)-min1); plot(V2(:,10)-min1);
%     legend({'hmean', 'amean', 'v1', 'v2'})
    
%     mr_column = (V1.*V2)./(V1+V2); %harmonic mean
%     mr_column = exp(log(V1) + log(V2) - log((V1+V2)) + log(2)); %less numerical error
else %phase 3
%     mrSingleColumn = repmat([1 1/3; 1/3 1], ceil(nSites_col/2), 1);
    % mrSingleColumn = repmat([1/3 1; 1 1/3], ceil(nSites_col/2), 1);    
    mrSingleColumn = repmat([1/3 1; 1 1/3].^.5, ceil(nSites_col/2), 1);    
%     mrSingleColumn = repmat([1 1/3; 1/3 1].^.5, ceil(nSites_col/2), 1);    
    mrSingleColumn = single(mrSingleColumn(1:nSites_col,:))';    
    mr_column = bsxfun(@times, single(mr(:, cviSite_track{1})), mrSingleColumn(1,:));
    mr_column = mr_column + bsxfun(@times, single(mr(:, cviSite_track{2})), mrSingleColumn(2,:));
end
end %func


function track_depth_(vcFile_prm)
% only for imec2 probe for now
% this won't work if internal reference used
if isempty(vcFile_prm)
    vcFile_prm = 'E:\Acute2\IMECIIiii2\10312014R13Block1128_all.prm'; 
end
[~, P] = file_prm2bin_(vcFile_prm);
[vcFile_bin, P] = imec_sites_(P);

% parameters
fMeanSubt = 0;

% input conditioning
nSamples_load = floor(P.sRateHz * P.tBin_track);  %load 10 sec bin
nBytes_load = nSamples_load * 2 * P.nChans;
nLoads = floor(file_bytes_(vcFile_bin) / nBytes_load);
if ~isempty(P.load_fraction_track)
    nLoads = round(nLoads * P.load_fraction_track);
end
vrTime_load = ((1:nLoads)-.5) * P.tBin_track;
sRateHz_lfp = P.sRateHz / P.nSkip_lfp;
[cviDist_depth, cvrPow_depth] = deal(cell(1, nLoads));
nSites = numel(P.viSite2Chan);
% nSites_valid = nSites - numel(P.viChan_ref);
eval(sprintf('hFun_track = @%s;', P.vcMode_track));
try eval(sprintf('clear %s;', P.vcMode_track)); catch; end

% fSingleColumn = 1; %experimental. interpolate a single line
if ~isempty(P.viDepth_track)
    P.viDepth_track(P.viDepth_track<1 | P.viDepth_track > floor(nSites/2)) = [];
end
fid = memmapfile(vcFile_bin, 'Format', {'int16', [P.nChans, nSamples_load], 'x'}, 'Offset', 0, 'Repeat', 1);
if ~P.fUseCache_track
    for iLoad=1:nLoads
        t1 = tic;        
        mrWav1 = (fid.Data.x); % load data
        if P.fUseLfp_track
            mrWav1 = mrWav1(:, 1:P.nSkip_lfp:end)'; %subsample and transpose
        else
            mrWav1 = mrWav1';
            sRateHz_lfp = P.sRateHz;
        end
%         mr1 = subsample_mean_(mr1, P.nSkip_lfp); % this averages the signal
        fid.Offset = fid.Offset + nBytes_load; %for the next loop   
        if ~isempty(P.viSites_dat)        
            mrWav1(:, P.viSites_dat) = mrWav1;    % dat 120 chan format
        else
            mrWav1 = mrWav1(:, P.viSite2Chan); %include reference channels
        end
    %     mr1(:, P.viChan_ref) = 0;
        mrWav1 = interp_ref_(mrWav1, P.viChan_ref);

        if fMeanSubt
            mrWav1 = single(mrWav1);
            mrWav1 = bsxfun(@rdivide, mrWav1, std(mrWav1)); %standardize
            mrWav1(:,P.viChan_ref) = 0;    
            mrWav1 = bsxfun(@minus, mrWav1, sum(mrWav1,2)/nSites_valid);     % mean subtract
            mrWav1(:,P.viChan_ref) = 0;
        end

        if P.fSingleColumn_track
            mr1_column = compute_column_(mrWav1, P.cviSite_track);
            [cviDist_depth{iLoad}, cvrPow_depth{iLoad}] = hFun_track(mr1_column, 'freqLim_excl', P.freqLim_excl_track, 'sRateHz', sRateHz_lfp, 'freqLim', P.freqLim_track);            
        else
            for iGroup=1:numel(P.cviSite_track)
                mr1_column = mrWav1(:, P.cviSite_track{iGroup});
                mrDist1 = hFun_track(mr1_column, 'freqLim_excl', P.freqLim_excl_track, 'sRateHz', sRateHz_lfp, 'freqLim', P.freqLim_track);
                if iGroup==1
                    trDist1 = zeros([size(mrDist1,1), size(mrDist1,2), numel(P.cviSite_track)]);                 
                end
                trDist1(:,:,iGroup) = mrDist1;
            end
            cviDist_depth{iLoad} = trDist1;
        end

        fprintf('Load %d/%d took %0.1fs\n', iLoad, nLoads, toc(t1));
    end %for
    save cviDist_depth cviDist_depth cvrPow_depth;
else
    load cviDist_depth;
end

vrShift = track_depth_shift_(cviDist_depth, P);

% plot
figure; plot(cell2mat(cvrPow_depth));
set(gcf,'Name', vcFile_bin, 'Color', 'w');
axis tight; grid on; xlabel('Electrode Depth'); ylabel('Power');

figure; hold on;
set(gcf,'Name', vcFile_bin, 'Color', 'w');
plot(vrTime_load, vrShift); 
if ~isempty(P.pix_per_sec_track) %descent experiment
    plot(vrTime_load([1 end]), [0, diff(vrTime_load([1 end]))] * P.pix_per_sec_track, 'r-');
end
axis tight; grid on; xlabel('Time (s)'); ylabel('# pixels shifted');
title(sprintf('Freq %0.1f~%0.1f Hz, maxSite_track:%d', P.freqLim_track, P.maxSite_track(end)), 'Interpreter', 'none');

S_shift = makeStruct(cviDist_depth, vrShift, P);
assignWorkspace(S_shift);

% display linear fit result
try
    [fitobject,gof] = fit((1:numel(S_shift.vrShift))', S_shift.vrShift, 'poly1');
    
    % estimate z resolution (intercept). 
%     n_per_step=10: del_z = 0.1400. smooth off, two-step off (95% CI)
%     n_per_step=5: del_z = 0.1448. smooth off, two-step off (95% CI)

    n_per_step = 10;
    n_steps = floor(numel(S_shift.vrShift)/n_per_step);
    mrShift = reshape(S_shift.vrShift(1:n_steps*n_per_step), n_per_step, []);
    mrShift = bsxfun(@minus, mrShift, mean(mrShift));
    if ~P.fSmooth_track, mrShift = smooth_121_mr_(mrShift); end
%     figure; plot(mrShift(:));
    mean(std(mrShift))
    std(std(mrShift))
%     std(mrShift(:))

%     mrShift_detrend = bsxfun(@minus, mrShift, S_shift.vrShift(1:n_per_step:end)');
%     [fitobject,gof] = fit((1:numel(mrShift_detrend))', mrShift_detrend(:), 'poly1')
catch
    disp(lasterr);
end
end %function


function vrDepth = track_depth_shift_(cviDist_depth, P)
vrDepth = imregcorr_matrix(cviDist_depth, P);

if P.fSmooth_track, vrDepth = smooth_121_mr_(vrDepth); end

end %func


function vrDepth = track_depth_shift_20160708(cviDist_depth, P)

vrShift1(iDepth) = imregcorr_matrix(cviDist_depth, P);

nDepths = numel(cviDist_depth);
mrShift = zeros(nDepths);
parfor iDepth=1:nDepths
    mrShift(:,iDepth) = track_depth_shift_pair_(...
        cviDist_depth{iDepth}, cviDist_depth, P); 
end
vrDepth = zeros(nDepths, 1);
vrShift_mid = mrShift(:, round(nDepths/2));
for iDepth = 1:nDepths
    vrDepth(iDepth) = median(vrShift_mid - mrShift(:,iDepth));
%     vrDepth(iDepth) = trimmean(vrShift_mid - mrShift(:,iDepth), 60);
end
vrDepth = vrDepth - vrDepth(1);

if P.fSmooth_track, vrDepth = smooth_121_mr_(vrDepth); end

end %end


function vrShift = track_depth_shift_old_(cviDist_depth, P)
% fSmooth_track = 1;
fTwoStep = 1;
% nw = 1;
% max_shift = 5;

%         [vrCorr1_odd, vnShift1] = xcorr_nan(cviDist_odd_depth{iDepth-1}, cviDist_odd_depth{iDepth}, 1);    
%         [vrCorr1_even, ~] = xcorr_nan(cviDist_even_depth{iDepth-1}, cviDist_even_depth{iDepth}, 1);    
%         vrCorr1 = vrCorr1_odd + vrCorr1_even;
%         [~, imax] = min(vrCorr1);
%         vrShift(iDepth) = vnShift1(imax);
         
% multi-depth interpolation
vrShift1 = track_depth_shift_pair_(cviDist_depth(1:end-1), cviDist_depth(2:end), P);
vrShift1 = cumsum(vrShift1);
if fTwoStep
    vrShift2 = track_depth_shift_pair_(cviDist_depth(1:end-2), cviDist_depth(3:end), P);
    vrShift2a = cumsum(vrShift2(1:2:end));
    vrShift2b = cumsum([vrShift1(1), vrShift2(2:2:end)]);

    % combine
    vrShift = zeros(numel(cviDist_depth), 1);
    vrShift(3:2:end) = (vrShift1(2:2:end) + vrShift2a)/2;
    vrShift(2:2:end) = (vrShift1(1:2:end) + vrShift2b)/2;
else
    vrShift = [0; vrShift1(:)];
end
if P.fSmooth_track, vrShift = smooth_121_mr_(vrShift); end

% figure; hold on; 
% vi1 = (1:numel(vrShift1))+1;
% vi2a = 3:2:numel(vrShift);
% vi2b = 2:2:numel(vrShift);
% plot(vi1, vrShift1, '.-', vi2a, vrShift2a, '.-', vi2b, vrShift2b, '.-', 1:numel(vrShift), vrShift, 'k.-');


% old, line based track
% if ~fPsd
%     cviDist_odd_depth{iLoad} = cellfun(@(mr)corr2line(mr, P.thresh_corr_track, P.viSites_track)', cviDist_odd_depth, 'UniformOutput', 0);
%     cviDist_even_depth{iLoad} = cellfun(@(mr)corr2line(mr, P.thresh_corr_track, P.viSites_track)', cviDist_even_depth, 'UniformOutput', 0);
%         
%     figure; set(gcf,'Name', vcFile_bin, 'Color', 'w');
%     subplot 211; imagesc(cell2mat(cviDist_odd_depth), [0 20]); xlabel('Depth (x 20um)'); ylabel('Dist (90% corr)'); axis xy; title('odd sites');
%     subplot 212; imagesc(cell2mat(cviDist_even_depth), [0 20]); xlabel('Depth (x 20um)'); ylabel('Dist (90% corr)'); axis xy; title('even sites');
% elseif ~fFft
%     figure; set(gcf,'Name', vcFile_bin, 'Color', 'w');
%     mr_odd = (cell2mat(cellfun(@(mr)mean(mr)', cviDist_odd_depth, 'UniformOutput', 0))); 
%     mr_even = (cell2mat(cellfun(@(mr)mean(mr)', cviDist_even_depth, 'UniformOutput', 0)));  
%     subplot 211; imagesc(mr_odd); xlabel('Session #'); ylabel('Site # (x 20 um)'); axis xy; title('odd sites');
%     subplot 212; imagesc(mr_even); xlabel('Session #'); ylabel('Site # (x 20 um)'); axis xy; title('even sites');
%     for iDepth = 1:numel(cviDist_odd_depth)
%         %@TODO: 2D xcorr
%         cviDist_odd_depth{iDepth} = mean(cviDist_odd_depth{iDepth})';
%         cviDist_even_depth{iDepth} = mean(cviDist_even_depth{iDepth})';
%     end
% else
    % fft cross covariance computation
%     vr1 = cviDist_even_depth{1};
%     vr2 = cviDist_even_depth{3};
%     figure; imagesc(abs(vr1(:,3:end)' * vr2(:,1:end-2)));  %max diag component after product
% end
end %function 


function vrShift1 = track_depth_shift_pair_(cviDist_depth1, cviDist_depth2, P)
vrShift1 = zeros(size(cviDist_depth2));
for iDepth = 1:numel(cviDist_depth2)  
    if iscell(cviDist_depth1)
        tr1 = cviDist_depth1{iDepth};
    else
        tr1 = cviDist_depth1;
    end
    tr2 = cviDist_depth2{iDepth};
    if ~isempty(P.viDepth_track)
        if ismatrix(tr1)
            tr1 = tr1(P.viDepth_track, P.viDepth_track); 
            tr2 = tr2(P.viDepth_track, P.viDepth_track); 
        else
            tr1 = tr1(P.viDepth_track, P.viDepth_track, :); 
            tr2 = tr2(P.viDepth_track, P.viDepth_track, :); 
        end
    end
    switch P.vcMode_track
        case 'mt_dpsd_mr'
            vrShift1(iDepth) = imregcorr_1D_xshift(tr1, tr2, P.thresh_corr_track, P.nw_lcm_track, P.max_shift_track);
        case {'mt_cpsd2_mr', 'mt_cpsd1_mr'}
%             vrShift1(iDepth) = imregcorr_1D_diag(tr1, tr2, P.maxSite_track, P.nw_lcm_track, P.max_shift_track);
            vrShift1(iDepth) = imregcorr_diag(tr1, tr2, P);
        otherwise
            vrShift1(iDepth) = imregcorr_1D(tr1, tr2, P.thresh_corr_track, P.nw_lcm_track, P.max_shift_track);
    end
end
end %func


function mr_mean = subsample_mean_(mr, nSkip)
% subsample columnn
% nChans = size(mr,1);
nSamples = floor(size(mr,2) / nSkip);
vi0 = (0:nSamples-1) * nSkip;

mr_mean = zeros(nSamples, size(mr,1), 'single');
for iSum = 1:nSkip
    mr_mean = mr_mean + single(mr(:, iSum + vi0)'); %subsample and transpose
end
mr_mean = mr_mean / nSkip; %average samples
end


function mrTf = imregcorr_thresh_(mrImg1, mrImg2, thresh)
mrImg1a = mr2triu(mrImg1);
mrImg2a = mr2triu(mrImg2);
mrImg1a(mrImg1a<thresh) = 0;
mrImg2a(mrImg2a<thresh) = 0;
mrImg1a(isnan(mrImg1a)) = 0;
mrImg2a(isnan(mrImg2a)) = 0;
tform = imregcorr(mrImg1a, mrImg2a, 'translation');
mrTf = tform.T;
end %func


function mr = smooth_121_mr_(mr)
mr = .25*mr([1,1:end-1],:) + .5*mr + .25*mr([2:end,end],:);
end %func


% function track_centroid_(vcFile)
% % plot all sites
% % After event detection
% [P, S0, Sclu] = get_param_(vcFile);
% iSite1 = Sclu.viSite_clu(iClu);
% viSite = P.miSites_pix(:, iSite1);
% viTime1 = Sclu.cviTime_clu{iClu};
% 
% [vrX_spk(viSpk1), vrY_spk(viSpk1), vrA_spk(viSpk1), mrVpp1_spk] = ...
%         spikePos_(viTime1, viSite1, S0.mrWav, P);
% 
% fPlotSite = 1;
% fPlotAllSites = 1;
% skip_spk = 1;
% fShowXpos = 0; fPlot_ampDist = 0;  
% fShiftPos = 1;  %show before and after difference
% 
% % viSites = [];
% [P, S0, Sclu] = get_param_();
% nSites = size(S0.mrWav,2);
% viClu_plot = [S0.iCluCopy, S0.iCluPaste];
% 
% % viSites = (1:30) + 30*(2-1); % comment out for local plot
% if fPlotSite %show all    
%     viClu_plot = []; %all clu    
%     if fPlotAllSites
%         viSites = 1:nSites;
%         skip_spk = 4;
%     else
%         iSite_clu = Sclu.viSite_clu(S0.iCluCopy); %only plot spikes from site     
%         viSites =  P.miSites_pix(:, iSite_clu)';
%     end    
% else
%     viSites = 1:nSites; 
%     skip_spk = 4;   fShowXpos = 0; fShiftPos = 1; fPlot_ampDist = 0;    
% end
% % vcFile_shift = 'D160709_descent_ph2o3.mat'; %correct shift, comment out to disable
% 
% % Compute   
% nSpikes = numel(Sclu.viTime);
% [vrX_spk, vrY_spk, vrA_spk, viClu_spk, vrT_spk] = deal(nan(nSpikes, 1, 'single'));
% % vl_spk = false(nSpikes, 1);
% cmrVpp_site = cell(size(viSites));
% for iSite = viSites
%     viSite1 = P.miSites_pix(:, iSite);
%     viSpk1 = find(Sclu.viSite == iSite);
%     viClu1 = Sclu.viClu(viSpk1);
%     if ~isempty(viClu_plot)        
%         vl1_clu = ismember(viClu1, viClu_plot);
%         viSpk1 = viSpk1(vl1_clu);
%         viClu1 = viClu1(vl1_clu);
%     end
%     if isempty(viSpk1), continue; end
%     viTime1 = Sclu.viTime(viSpk1);    
%     [vrX_spk(viSpk1), vrY_spk(viSpk1), vrA_spk(viSpk1), mrVpp1_spk] = ...
%         spikePos_(viTime1, viSite1, S0.mrWav, P);
%     viClu_spk(viSpk1) = viClu1;
%     vrT_spk(viSpk1) = viTime1;
% %     vl_spk(viSpk1) = 1;
% %     mrVpp2_spk = zeros(nSites, size(mrVpp1_spk,2), 'single');
% %     mrVpp2_spk(viSite1, :) = mrVpp1_spk;
%     if fPlot_ampDist, cmrVpp_site{iSite} = sort(mrVpp1_spk, 'descend'); end
% end %for
% 
% % select
% vl_spk = ~isnan(vrY_spk);
% vrA_spk = (vrA_spk(vl_spk));
% vrY_spk = vrY_spk(vl_spk) / P.um_per_pix;
% vrX_spk = vrX_spk(vl_spk) / P.um_per_pix;
% vrT_spk = single(vrT_spk(vl_spk)) / P.sRateHz;
% viClu_spk = viClu_spk(vl_spk);
% if fPlot_ampDist, plot_ampDist_(cmrVpp_site, P); end
% 
% % centroid spike plot time bin
% if fShiftPos
%     [vrT_shift, vrY_spk, vrY_shift] = shift_spike_(vrT_spk, vrY_spk, vrA_spk, P);
% end
% 
% % vrA_spk = log10(vrA_spk); %plot purpose
% if fPlotSite
%     % sort by ampl and keep half
%     [vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk] = ...
%         sort_ascend_(vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk);
% %     viSelect = round(numel(vrA_spk)/2):skip_spk:numel(vrA_spk); %top half subsample
%     viSelect = 1:skip_spk:numel(vrA_spk); %top half subsample
%     [vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk] = ...
%         select_vr_(vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk, viSelect);    
% %     vrY_spk = shift_spk_pos_(vcFile_shift, vrT_spk, vrY_spk);
%     ylim_pix = [floor(min(viSites / 2)), ceil(max(viSites / 2))];
% else
%     ylim_pix = round(median(vrY_spk)) + [-1,1] * floor(P.maxSite);
% end
% 
% %------------
% % Draw
% hFig = resize_figure([], [.5 0 .5 .5]); clf; %2001
% try set(hFig, 'Name', P.vcFile_prm); catch; end
% AX = [];
% if isempty(S0.iCluPaste) || fPlotSite % Show one unit or one site
%     if fShowXpos
%         AX(1) = subplot(211); hold on;
%     else
%         AX = gca; hold on;
%     end
%     set(gca, 'YLim', ylim_pix);
%     copy_axes_x_(gca, P.hFigTime);
%     colormap(flipud(colormap('gray')));     
%     ylabel('Y pos [pix]');
%     vcTitle = 'Color: log10 Vpp [uV]';
%     if fPlotSite
%         title(sprintf('%s; Site%d', vcTitle, viSites(1)));
%     else
%         title(sprintf('%s; Clu%d', vcTitle, S0.iCluCopy));
%     end
%     colorbar(gca);
%     scatter(vrT_spk, vrY_spk, 5, log10(vrA_spk), 'filled'); grid on;
%     
%     if ~fShowXpos, return; end
%     AX(2) = subplot(212); hold on;
%     set(gca, 'YLim', [0 1.5]);
%     copy_axes_x_(gca, P.hFigTime);
%     colormap(flipud(colormap('gray')));     
%     xlabel('Time (s)'); ylabel('X pos [pix]');
%     colorbar(gca);   
%     set(AX, 'CLim', [.5 3.5]);
%     scatter(vrT_spk, vrX_spk, 5, log10(vrA_spk), 'filled'); grid on;
% else %show black and red
%     if fShowXpos
%         AX(1) = subplot(211); hold on;
%     else
%         AX = gca; hold on;
%     end
%     set(gca, 'YLim', ylim_pix);
%     copy_axes_x_(gca, P.hFigTime);
%     colormap(gca, [1 0 0; 0 0 0]); 
%     colorbar(gca, 'off');    
%     ylabel('Y pos [pix]');
%     title(sprintf('Clu%d: black; Clu%d: red', S0.iCluCopy, S0.iCluPaste));
%     scatter(vrT_spk, vrY_spk, 5, viClu_spk, 'filled'); grid on;        
%     
%     if ~fShowXpos, return; end
%     AX(2) = subplot(212); hold on;
%     set(gca, 'YLim', [0 1.5]);
%     copy_axes_x_(gca, P.hFigTime);
%     colormap(gca, [1 0 0; 0 0 0]); 
%     colorbar(gca, 'off');
%     xlabel('Time (s)'); ylabel('X pos [pix]');
%     scatter(vrT_spk, vrX_spk, 5, viClu_spk, 'filled'); grid on;           
% end
% % linkaxes(AX, 'x');
% % mouse_figure(hFig);
% end %func


function download_sample_(vcMode)
S_cfg = read_cfg_();
if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
switch lower(vcMode)   
    case {'sample', 'neuropix2', 'neuropixels2', 'phase2', 'phaseii'}
        csLink = S_cfg.path_sample_phase2;
    case {'sample3', 'neuropix3' 'neuropixels3', 'phase3', 'phaseiii'}
        csLink = S_cfg.path_sample_phase3;
    otherwise
        disp('Invalid selection. Try "jrclust download sample".');
        return;
end %switch

t1 = tic;
fprintf('Downloading sample files. This can take up to several minutes.\n');
download_files_(csLink);
fprintf('\tAll sample files downloaded. Took %0.1fs\n', toc(t1));
end


function download_files_(csLink, csDest)
% download file from the web
if nargin<2, csDest = link2file_(csLink); end

for i=1:numel(csLink)    
    % download from list of files
    fprintf('\tDownloading %s: ', csLink{i});
    vcFile_out1 = websave(csDest{i}, csLink{i});
	fprintf('saved to %s\n', vcFile_out1);
end %for
end


function csFile = link2file_(csLink)
csFile = cell(size(csLink));
for i=1:numel(csLink)        
    vcFile1 = csLink{i};
    iBegin = find(vcFile1=='/', 1, 'last'); % strip ?    
    if ~isempty(iBegin), vcFile1 = vcFile1(iBegin+1:end); end

    iEnd = find(vcFile1=='?', 1, 'last'); % strip ?
    if ~isempty(iEnd), vcFile1 = vcFile1(1:iEnd-1); end
    csFile{i} = vcFile1;
end
end %func


function [P, S0, Sclu] = get_param_(vcFile_prm)
if nargin>0
    P = loadParam_(vcFile_prm);
    Sclu = load_clu_(P.vcFile_prm);
    mrWav = loadWavFile(P.vcFile_prm);    
    S0 = makeStruct(P, Sclu, mrWav); 
    set(0, 'UserData', S0);
else
    % @TODO: use this over workspace-based parameter caching
    S0 = get(0, 'UserData');    
    P = S0.P; %get parameter from user data
    Sclu = S0.Sclu;
end
end %func


function S = load_(vcFile_prm, vcType)
if nargin<2, vcType = []; end
S0 = get(0, 'UserData');
if isempty(vcType)
    S = S0;     
else
    try
        S = S0.(vcType);
    catch
        S = [];
    end
end
% switch lower(vcType)
%     case 'sevt'
%     case 'sclu'
% 
% end
end %func


function val = read_cfg_(vcName)
% read configuration file that stores path to folder
% @TODo: cross-platform support can be fixed here
% load from default.cfg but override with user.cfg if it exists
if ~exist('default.cfg', 'file')
   S_cfg = struct( ...
        'path_dropbox', 'C:\Dropbox (HHMI)\Git\jrclust', ...
        'path_backup', 'C:\backup', ...
        'path_alpha', 'C:\Dropbox (HHMI)\Git\jrclust_alpha', ...
        'sync_list', {'default.prm', '*.txt', '*.m', '*.ptx', '*.cu', '*.prb', 'default.cfg', 'JRClust manual.docx', ...
        'path_sample_phase2', 'https://www.dropbox.com/s/t0my3nf8dmpzmr6/sample.meta?dl=1', 'https://www.dropbox.com/s/xhgbb624h2tls6h/sample.bin?dl=1', ... 
        'path_sample_phase3', 'https://www.dropbox.com/s/9vuldabqw68ilpa/sample3.meta?dl=1', 'https://www.dropbox.com/s/lcowg4c9f4xiat5/sample3.bin?dl=1', ...
        });
else
    S_cfg = file2struct('default.cfg');
end
try
    S_cfg1 = file2struct('user.cfg');
    S_cfg.path_dropbox = S_cfg1.path_dropbox;
    S_cfg.path_backup = S_cfg1.path_backup;
catch
    disperr();
end
if nargin==0
    val = S_cfg; 
else
    val = S_cfg.(vcName);
end
end %func


% function [vr_sum, vr_diff] = sum_diff_(vr1, vr2)
% vr_sum = (vr1+vr2);
% vr_diff = (vr1-vr2);
% end %func


function [fConst, fOdd, fEven] = cable_test_phase3_(mnWav)
% vrMax = double(max(mnWav));
% vrMin = double(min(mnWav));
% vrSum = vrMax + vrMin;
% vrDiff = vrMax - vrMin;

[vlWav_time, vlWav_chan] = cable_test_phase3_match_(mnWav);
fConst = all(vlWav_time); 
fOdd = all(vlWav_chan(1:2:end));
fEven = all(vlWav_chan(2:2:end));

% fOdd = std(vrSum(1:2:end)) == 0;
% fEven = std(vrSum(2:2:end)) == 0;
end


function [vlWav_time, vlWav_chan] = cable_test_phase3_match_(mnWav)
mlWav = true(size(mnWav));
mlWav(:, 1:2:end) = mnWav(:, 1:2:end) == -170;
mlWav(:, 2:2:end) = mnWav(:, 2:2:end) == 170;
vlWav_chan = all(mlWav);
vlWav_time = all(mlWav,2);
end


function cable_test_(vcFile)
% cable test for phase 3
% all odd must be the same value, even must be the same, no change over time
% wild card
fSortByDate = 1; %sort files by date created

csFiles = dir_file_(vcFile, fSortByDate);
P = struct('viRef', [], 'nSamples_load', 20e6); %load 15.4GB approx.
for iFile = 1:numel(csFiles)
%     vcFile_bin = sprintf('%s\\20160715-1630_CableTest_g0_t%d.imec.ap.bin', vcDir, iFile);
% %     fprintf('loading: %s: ', vcFile_bin);
    vcFile_bin = csFiles{iFile};
    [mnWav, mnLfp, mnAux, Smeta] = import_imec3(vcFile_bin, P);    
    [fConst_ap, fOdd_ap, fEven_ap] = cable_test_phase3_(mnWav);
    [fConst_lf, fOdd_lf, fEven_lf] = cable_test_phase3_(mnLfp);
    vlTest = [fConst_ap, fOdd_ap, fEven_ap, fConst_lf, fOdd_lf, fEven_lf];
    fSuccess = all(vlTest);
    if fSuccess
        fprintf('\t%s: PASS.\n', vcFile_bin);
    else
        fprintf('\tFAIL: fConst_ap:%d, fOdd_ap:%d, fEven_ap:%d, fConst_lf:%d, fOdd_lf:%d, fEven_lf:%d\n', int32(vlTest));
%         [vlWav_ap_time, vlWav_ap_chan] = cable_test_phase3_match_(mnWav);
%         [vlWav_lf_time, vlWav_lf_chan] = cable_test_phase3_match_(mnLfp);
%         fprintf('\tAP-time:%f; AP-chan:%f; LF-time:%f; LF-chan:%f\n', ...
%             mean(vlWav_ap_time), mean(vlWav_ap_chan), mean(vlWav_lf_time), mean(vlWav_lf_chan));
    end
end
assignWorkspace(mnWav, mnLfp, mnAux, Smeta);
end %func


function fSuccess = compile_cuda_(S_cfg)
if nargin<1, S_cfg = read_cfg_(); end
t1 = tic;
csFiles_cu = S_cfg.csFiles_cu;
disp('Compiling CUDA codes...');
fSuccess = 1;
for i=1:numel(csFiles_cu)
    vcCmd1 = sprintf('nvcc -ptx %s', csFiles_cu{i});
    fprintf('\t%s\n', vcCmd1);
    status = system(vcCmd1);
    fSuccess = fSuccess && (status==0);
end
fprintf('\tFinished compiling, took %0.1fs\n', toc(t1));
end %func


function install_jrclust_()
% compile nvcc codes and create a user.cfg
try
    S_cfg = read_cfg_();
    update_version_(S_cfg);
    if ~install_cuda_(), error('CUDA installation error'); end
    if ~install_vc_(), error('Visual Studio installation error'); end
    if ~compile_cuda_(S_cfg), error('CUDA compile error'); end
    copyfile('template.cfg', 'user.cfg', 'f');
    edit('user.cfg');
    msgbox('Edit the user''s configuration if needed.', 'modal');
catch
    disp(lasterr());
    fprintf(2, 'Installation failed. Restart Matlab and run installation again.\n'); 
end
end %func


function fSuccess = install_vc_()
% check for visual c++ installation
fSuccess = 0;
vcLink_vc = 'https://www.microsoft.com/en-us/download/details.aspx?id=34673';

if exist('C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin', 'dir')
    add_path_('C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\bin');
    fSuccess = 1;
elseif exist('C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin', 'dir')
    add_path_('C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin');
    fSuccess = 1;
else
    fprintf(2, 'Visual studio is not installed. Please download and install Visual Studio Express 2012 from the link below:\n');
    disp(['<a href = "', vcLink_vc, '">Visual studio download.</a>']);
end
end %func


function fSuccess = install_cuda_()
% Sets the system path
fSuccess = 0;

% check for CUDA installation
try
    S_gpu = gpuDevice();
    if isempty(S_gpu), error('CUDA-compatible GPU is not found.'); end
catch
    fprintf(2, 'CUDA-compatible GPU is not found. Recommended GPU: GTX 980 Ti 6GB or Titan X 12 GB\n');
    return;
end

% Set CUDA path
vcPath_cuda = ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v', sprintf('%0.1f', S_gpu.ToolkitVersion)];
if ~exist(vcPath_cuda, 'dir')
    fprintf(2, 'CUDA Toolkit is not installed. Download from the link below.\n');
    disp(['<a href = "', getCudaLink_(S_gpu), '">CUDA Toolkit download.</a>']);
    return;
end
csPath_cuda = cellfun(@(vc)[vcPath_cuda, filesep(), vc], {'libnvvp', 'bin'}, 'UniformOutput', 0);
add_path_(csPath_cuda);

% set GPU watchdog timer value
val = get_registry_('HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers', 'TdrDelay');
if isempty(val) || val ~= 0
    set_registry_('HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers', 'TdrDelay', 0);
end

% set CUDA cache mazsize
% https://github.com/FORTH-ModelBasedTracker/HandTracker/issues/7
set_env_('CUDA_CACHE_MAXSIZE', 1073741824);
set_env_('CUDA_CACHE_DISABLE', 0);
% get_env_('CUDA_CACHE_MAXSIZE')
fSuccess = 1;
end %func


function val = get_registry_(vcKey, vcName)
% vcKey = 'HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers';
% vcName = 'TdrDelay';

if ~ispc(), fprintf(2, 'Registry only works for Windows.'); val = []; return; end

iEnd = find(vcKey=='\', 1, 'first');
if isempty(iEnd), val=[]; return; end
vcRootKey = vcKey(1:iEnd-1);
vcSubKey = vcKey(iEnd+1:end);
try
    val = winqueryreg(vcRootKey, vcSubKey, vcName);
catch
    val = [];
end
end %func


function fSuccess = set_registry_(vcKey, vcName, vcValue)
fSuccess = 0;

if ~ispc(), fprintf(2, 'Registry setting only works for Windows.'); return; end
% vcKey = 'HKEY_LOCAL_MACHINE\System\CurrentControlSet\Control\GraphicsDrivers';
% vcName = 'TdrDelay';
% vcValue = '0';

regFileName = 'c:\\jrclust_temp_reg.reg';
fp = fopen(regFileName,'wt');
% fprintf(fp,'REGEDIT4\n');
fprintf(fp, 'Windows Registry Editor Version 5.00\n\n');
fprintf(fp,'[%s]\n',vcKey);
if ischar(vcValue)
    vcValue = strrep(vcValue,'\','\\');     %  escape backslashes
    fprintf(fp,'"%s"="%s"\n', vcName, vcValue);
%     fprintf(fp,'%s%s%s%s%s%s\n','"', vcName, '"=', '"', vcValue, '"' );
    fSuccess = 1;
else
    if numel(vcValue)==1
        fprintf(fp,'"%s"=dword:%.8X', vcName, int32(vcValue)); % write as hexidecimal
    end
    fSuccess = 1;
end
fclose(fp);

if fSuccess    
    fSuccess = system(sprintf('C:\\Windows\\regedit.exe /S "%s"', regFileName));
%     open(regFileName);
    fSuccess = fSuccess == 0;
end

drawnow;
delete(regFileName);
end %func


function fSuccess = set_env_(valname, value)
% fSuccess = set_registry_('HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment', valname, value);
if ischar(value)
    fSuccess = system(sprintf('setx %s %s', valname, value));
else
    fSuccess = system(sprintf('setx %s %d', valname, int32(value)));
end
end %func


function val = get_env_(valname)
% val = get_registry_('HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment', valname);
val = getenv(valname);
end %func


function add_path_(csPath_cuda)
if ischar(csPath_cuda), 
    csPath_cuda={csPath_cuda}; 
end

% add following path to the system, in front of the path
vcPath = getenv('PATH');
pathLen = numel(vcPath);
for i=1:numel(csPath_cuda)
    if isempty(strfind(vcPath, csPath_cuda{i}))
        vcPath = [csPath_cuda{i}, ';', vcPath];
    end
end
if numel(vcPath) > pathLen
    % setenv('PATH', vcPath); %the change does not persist
%     set_registry_('HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\Session Manager\Environment', 'PATH', vcPath);
    set_env_('PATH', vcPath);
end
end %func


function vcLink = getCudaLink_(S_gpu)
if nargin<1, S_gpu = gpuDevice(); end
if S_gpu.ToolkitVersion < 7.5    
    vcLink = sprintf('https://developer.nvidia.com/cuda-toolkit-%d', round(S_gpu.ToolkitVersion*10));
else
    vcLink = 'https://developer.nvidia.com/cuda-downloads';
end
end %func


function [csFiles, viValid] = filter_files_(csFiles)
% Files must exist and non-zero byte
if ischar(csFiles)
    if any(csFiles=='*')
        csFiles = dir_file_(csFiles, 1); %sort by dates
    else
        csFiles = {csFiles}; %make it a cell
    end
end

% filter based on file presence and bytes
vlValid = false(size(csFiles));
for iFile=1:numel(csFiles)    
    S_dir1 = dir(csFiles{iFile});
    if isempty(S_dir1), continue; end
    if S_dir1.bytes == 0, continue; end
    vlValid(iFile) = 1;
end
viValid = find(vlValid);
csFiles = csFiles(viValid);
end %func


function track_set_(vcFile_set)
% Analyze a collection of clustered files
% Based on S160615_batch_imec3.m

if ~matchFileExt(vcFile_set, '.set'), msgbox('Must provide mycollection.set file'); return ;end
edit(vcFile_set); %show set file currently working

% Load S.set file
S_set = file2struct('default.set');
S_set = appendStruct(S_set, file2struct(vcFile_set));

% Filter files
csFiles_set = filter_files_(S_set.csFiles_set); %keep valid files only. 
csFiles_set = excludeString(csFiles_set, S_set.csFile_excl);
% @TODO: include a date filtering and name filtering . combine in filter_files
if isempty(csFiles_set), msgbox('No files found', 'modal'); end
msgbox(csFiles_set, 'modal');

[mrAmp90_site, mrRate_site, mrSnr_site, mrThresh_site] = deal(zeros(numel(csFiles_set), 3)); 
vnClu_file = zeros(numel(csFiles_set), 1);
vrDate_file = nan(size(csFiles_set));

% load directory
for i1=1:numel(csFiles_set)
    vcFile_prm1 = csFiles_set{i1};
    try        
        S_prm1 = load_session(vcFile_prm1, S_set);
        if isempty(S_prm1), continue; end
        [mrAmp90_site(i1,:), mrRate_site(i1,:), mrSnr_site(i1,:), mrThresh_site(i1,:)] =  ...
            multifun_(@(x)quantile(x, [.25,.5,.75]), ...
                S_prm1.vrAmp90_site, S_prm1.vrRate_site, S_prm1.vrSnr_site, S_prm1.vrThresh_site);
        vrDate_file(i1) = S_prm1.date_file;
        vnClu_file(i1) = S_prm1.nClu;
        disp(vcFile_prm1);
    catch
        fprintf('%s: %s\n', vcFile_prm1 , lasterr());
    end        
end %iFile1
viKeep = find(~isnan(vrDate_file));
[mrAmp90_site, mrRate_site, mrSnr_site, vrDate_file, vnClu_file, mrThresh_site] = ...
    select_vr_(mrAmp90_site, mrRate_site, mrSnr_site, vrDate_file, vnClu_file, mrThresh_site, ...
        viKeep);
if isempty(S_set.vcStart_date)
    vrDate_postop = vrDate_file(:) - vrDate_file(1);
else
    vrDate_postop = vrDate_file(:) - datenum(S_set.vcStart_date);
end
[vcDir1, ~, ~] = fileparts(csFiles_set{1});
vcTitle = sprintf('%s', vcDir1);

%----------------------
% Plot
figure; ax=[];  set(gcf, 'Name', vcFile_set, 'Color', 'w');

ax(end+1)=subplot(511); plot(vrDate_postop, mrAmp90_site, 'o-'); 
ylabel('Amp (uV, 90th)'); grid on; ylim([0 200]);
title(vcTitle, 'Interpreter', 'none');

ax(end+1)=subplot(512); plot(vrDate_postop, mrRate_site, 'o-'); 
ylabel('Rate (Hz)');  grid on; ylim([.1 1e2]); set(gca, 'YScale', 'log');

ax(end+1)=subplot(513); plot(vrDate_postop, vnClu_file, 'o-'); 
ylabel('# Clusters');  grid on; 

ax(end+1)=subplot(514); plot(vrDate_postop, mrThresh_site, 'o-'); 
ylabel('Threshold (uV)');  grid on; 

ax(end+1)=subplot(515); plot(vrDate_postop, mrSnr_site, 'o-'); 
ylabel('SNR');  grid on; ylim([0 20]); ylim([0 15]);

xlabel('Days since implantation'); 
linkaxes(ax, 'x');
end %func


function varargout = select_vr_(varargin)
% sort ascend
viKeep = varargin{end};
if islogical(viKeep), viKeep = find(viKeep); end
for i=1:(nargin-1)
    if isvector(varargin{i})
        varargout{i} = varargin{i}(viKeep);
    else
        varargout{i} = varargin{i}(viKeep, :);
    end
end
end %func


function varargout = multifun_(hFun, varargin)
if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    varargout{i} = hFun(varargin{i});
end
end


function Sevt = getVpp_Sevt_(Sevt, mrWav, P)
% get feature from here
% Calculate Vpp per site where spike centers occur
% [mrVpp, mrEnergy, viSite] = getVpp_Sevt_(mrWav, viSpk, viSite, P)
% later modify viSite based on the centroid

fprintf('Calculating Vpp '); t1=tic;
nChans = size(mrWav, 2);
% [cviSpk_site, cmrVpp_site] = deal(cell(nChans, 1));
Sevt.cviSpk_site = cell(nChans, 1);
miSite_vpp = findNearSites(P.mrSiteXY, P.maxSite * 2); % add a wiggle room?
nSites_vpp = size(miSite_vpp, 1);
nSpk = numel(Sevt.viSpk);
Sevt.mrVpp_spk = zeros(nSites_vpp, nSpk, 'single'); % vpp feature
% Sevt.tmrPc_spk = zeros(P.nPcPerChan, 
% get three principal vector per each
for iChan = 1:nChans %try parfor?
   viSpk1 = find(Sevt.viSite_spk == iChan); %find spikes centered at certain channel
   if isempty(viSpk1), continue; end
   trWav1 = mr2tr2(mrWav, Sevt.viSpk(viSpk1), P.spkLim, 0, miSite_vpp(:, iChan));      
   Sevt.mrVpp_spk(:, viSpk1) = single(max(trWav1) - min(trWav1)) * P.uV_per_bit;
%    cmrEnergy_site{iChan} = std(single(trWav1)) * P.uV_per_bit;
   fprintf('.');
end
fprintf(' took %0.1fs\n', toc(t1)); 

% output
% Sevt.cviSpk_site = cviSpk_site;
% Sevt.cmrVpp_site = cmrVpp_site;
% Sevt.miSite_vpp = miSite_vpp;
% Sevt.miSite_vpp = miSite;
end %func


function Sevt = centroid_Sevt_(Sevt, mrWav, P)
% determine spike centroid and shift the center

nSites = numel(P.viSite2Chan);
miSite = findNearSites(P.mrSiteXY, P.maxSite);
Sevt.mrPosXY_spk = zeros(numel(Sevt.viSpk), 2, 'single');
Sevt.viSite_spk = zeros(size(Sevt.viSite), 'like', Sevt.viSite); %nearest site
for iSite = 1:nSites
    % Todo: this could increase
    viSite1 = miSite(:, iSite);
    mrSiteXY1 = P.mrSiteXY(viSite1, :);
    viSpk1 = find(Sevt.viSite == iSite);
    mrVpp1 = mr2tr2(mrWav, Sevt.viSpk(viSpk1), P.spkLim, 0, viSite1);
    mrVpp1 = single(squeeze(max(mrVpp1) - min(mrVpp1))) * P.uV_per_bit;
    mrPosXY1 = centroid_mr_(mrVpp1, mrSiteXY1);
    Sevt.mrPosXY_spk(viSpk1, :) = mrPosXY1;
    [~, viSite_nearest1] = min(eucldist2_(mrSiteXY1, mrPosXY1));
    Sevt.viSite_spk(viSpk1) = viSite1(viSite_nearest1);
end
end %func


function [mrWav_spk, viSpk1] = Sevt_wav_site_(mrWav, Sevt, iSite, P)
% return spike waveform at certain site
% subtract mean
if nargout>=2
    viSpk1 = find(Sevt.viSite == iSite);
    viTime1 = Sevt.viSpk(viSpk1);
else
    viTime1 = Sevt.viSpk(Sevt.viSite == iSite);
end
% mean subtract and single
mrWav_spk = mr2tr3(mrWav, P.spkLim, viTime1, iSite, 0); %no reshape needed
% mrWav_spk = squeeze(mr2tr2(mrWav, viTime1, P.spkLim, 0, iSite));
% mrWav_spk = single(mrWav_spk);
end %func


function [vrPv, mrPv, vrPv_mean, mrWav_spk0, Sevt_subset] = pca_detect_(mrWav, P)
% determine spike centroid and shift the center
% find a matched filter to run deconvolution
fRedo_pca = 1;
nPc_max = 6;
fPcaMode = 4;  % 3: kmeans, 1 for pca; 4: kurtosis
fMatchMode = 2; %2 for mean, 1 for pca
fWhiten = 0;
fMean_zscore = 1;
fAlignMean = 1;

t1 = tic;
fprintf('Determining matched filter kernel\n');
TEMPLATE_DURATION = 10; %use 10 sec template duration by default

% Do a preliminary detection to determine pca
P.nlim_detect = [1, round(P.sRateHz * TEMPLATE_DURATION)];
P.nlim_detect(2) = min(size(mrWav,1), P.nlim_detect(2));

P.vcSpatialFilter = 'none'; %no need to do spatial sum in this case
P.spkThresh_uV = [];

% Collect spike waveforms
for i=1:(fRedo_pca + 1)
    if i==1, P.vrPv = []; end %clear matched filter
    Sevt_subset = spikeDetect_Sevt_(mrWav, P);
    nSites = numel(P.viSite2Chan);
    cmrPv = cell(1, nSites);
    mrWav_spk = cell(1, nSites);
    for iSite = 1:nSites
        mrWav_spk1 = Sevt_wav_site_(mrWav, Sevt_subset, iSite, P);   
%         mrWav_spk{iSite} = mrWav_spk1;
        if isempty(mrWav_spk1), continue; end
        try
            mrWav_spk1 = single(mrWav_spk1);
%             mrWav_spk1 = zscore(mrWav_spk1);
%             mrWav_spk1 = bsxfun(@minus, mrWav_spk1, mean(mrWav_spk1));
            if fPcaMode==2
                [~, cmrPv{iSite}] = pca(mrWav_spk1, 'NumComponents', nPc_max);  
            end
            mrWav_spk{iSite} = mrWav_spk1;
        catch
            ;
        end
    end
%     mrWav_spk = single(cell2mat_(mrWav_spk));
    mrWav_spk0 = cell2mat_(mrWav_spk);
    mrWav_spk = mrWav_spk0;
    if fMean_zscore
        mrWav_spk = zscore(mrWav_spk);
    end
    if fAlignMean
        mrWav_spk = align_mean(mrWav_spk, -2:2, 3);
    end
    vrPv_mean = mean(mrWav_spk, 2);
    vrLat = [];
    if fPcaMode == 1
        % matlab way
%         [mrPv0, mrPc0, vrLat0] = pca(mrWav_spk', 'NumComponents', nPc_max, 'Algorithm', 'eig');
%         vrLat0 = vrLat0(1:nPc_max) / sum(vrLat0(1:nPc_max));
        
        % my way
        mrWav_spk = bsxfun(@minus, mrWav_spk, mean(mrWav_spk)); % mean subtract
        [mrPv, mrPc1, vrLat, vrPv_mean1] = pca_eig(mrWav_spk);
%         mrPc2 = zscore(mrPv(:,1:nPc_max))' * bsxfun(@minus, mrWav_spk, vrPv_mean);
        mrPv = zscore(mrPv(:, 1:nPc_max)); 
    elseif fPcaMode == 4
        mrPv = cluster_kurtosis(mrWav_spk, P.nPcPerChan);
    elseif fPcaMode == 3 %k means
        if P.nPcPerChan == 1
            mrPv = vrPv_mean;
        else
            viKmeans = kmeans(mrWav_spk', P.nPcPerChan);
            mrPv = zeros(size(mrWav_spk,1), P.nPcPerChan, 'single');
            for iPc=1:P.nPcPerChan
                mrPv(:,iPc) = median(mrWav_spk(:, viKmeans == iPc), 2);
            end            
        end        
        mrPv = zscore(mrPv);
    else
    %     mrWav_spk = norm_mr(mrWav_spk);
    %     mrPvPv = cell2mat_(cmrPv);
        trPvPv = cell2mat_(cmrPv);
        trPvPv = reshape(trPvPv, size(trPvPv,1), nPc_max, []);
        mrPvPv = squeeze(trPvPv(:,1,:));
        vrPv_mean = mean(mrPvPv, 2);
    %     [~, ~, vrLat] = mrPvPv_to_mrPv_(mrPvPv, P);
        [~, ~, vrLat] = pca(mrPvPv);
        mrPv = cell(1, nPc_max);
        for iPc=1:nPc_max
            mrPv{iPc} = pca(squeeze(trPvPv(:,iPc,:))', 'NumComponents', 1);
        end    
        mrPv = cell2mat(mrPv);
    end
%     mrPv = bsxfun(@minus, mrPv, mean(mrPv));
%     mrPv = norm_mr(mrPv);    
%     mrPv = zscore(mrPv);
%     mrPv = mrPv * diag(1./sqrt(vrLat)); %whiten
   
%     [mrPv, vrLat, vrPv_mean] = pca_jjj(mrWav_spk, 'nPc_pca', P.nPcPerChan, 'nSamples_pca', 10000, 'fMeanSubt_pca', 1);
    if ~isempty(vrLat)
        vrLat = vrLat(1:nPc_max)' / sum(vrLat(1:nPc_max));
        if fWhiten
            mrPv = mrPv * diag(1./sqrt(vrLat));
        end
    end
%     mrWav_spk = subsample_mr(mrWav_spk, 10000);
%     vrPv_mean = mean(mrWav_spk,2);
%     mrWav_spk = bsxfun(@minus, mrWav_spk, vrPv_mean);
%     [mrPv, vrLat] = eig(mrWav_spk*mrWav_spk');
%     vrLat = flipud(diag(vrLat));
%     mrPv = mrPv(:, end:-1:end-P.nPcPerChan+1);    
%     mrPc = mrPv' * mrWav_spk;
    viPcFlip = find(mrPv(1-P.spkLim(1),:) >0);
    mrPv(:,viPcFlip) = -mrPv(:,viPcFlip);
        
    % convert to projection matrix
%     mrPv = pinv(mrPv');
%     mrPv = mrPv(:, 1:nPc_max);


%     if ~isempty(mrWav_spk2) %flip polarity so that pc is positive
%         viFlip = find(mean(mrWav_spk2' * mrPv) < 0);
%         mrPv(:,viFlip) = -mrPv(:,viFlip);
%     end
    if fMatchMode==1
        % matched filter
        vrPv = flipud(-zscore(mrPv(:,1)));
    else
        vrPv = flipud(-zscore(vrPv_mean));
    end
    P.vrPv = vrPv;
    if ~isempty(vrLat)
        fprintf('\n\tPC1 explains %0.1f%% of variance; vrLat:\n', vrLat(1)*100);
        fprintf('\n\tlat(%%): %s\n', sprintf('%0.1f ', cumsum(vrLat)*100));
    end
end
fprintf('\n\tDetermining matched filter took %0.1fs.\n', toc(t1));
end %func


function [mrPv, vrPv, vrLat] = mrPvPv_to_mrPv_(mrPvPv, P)
fWhiten = 0;

imin0 = 1-P.spkLim(1);
vrPv_min = mrPvPv(imin0,:);
viFix = find(vrPv_min >= 0);
vrPv_min(viFix) = min(mrPvPv(:,viFix));
mrPvPv = bsxfun(@rdivide, mrPvPv, -vrPv_min); %normalize by min
if fWhiten    
    [mrPv, vrLat] = zca(mrPvPv, 'fPca_zca', 1, 'nPc_zca', P.nPcPerChan);
else
    [~, mrPv, vrLat] = pca(mrPvPv, 'NumComponents', P.nPcPerChan);
end
mrPv = single(mrPv);
for iPc = 1:P.nPcPerChan
    if mrPv(imin0, iPc)>0
        mrPv(:,iPc) = -mrPv(:,iPc);
    end
end
mrPv = bsxfun(@minus, mrPv, mean(mrPv)); %not necessary
vrPv = -flipud(mrPv(:,1));
if ~fWhiten
    mrPv = bsxfun(@rdivide, mrPv, sqrt(sum(mrPv.^2))); %do not scale
end
vrPv = vrPv / sqrt(sum(vrPv.^2)); %no need
% [mrPv, vrLat] = zca(mrPvPv, 'fPca_zca', 0, 'nPc_zca', P.nPcPerChan);
end %func

% mrWav1 = zeros(size(mrWav), 'single');
% for iSite = 1:nSites
%     mrWav1(:,iSite) = matched_filt_(mrWav(:,iSite), vrPv, P.spkLim);
% end
% 
% mrWav = single(mrWav); mrWav = mrWav / std(mrWav(:));
% mrWav1 = single(mrWav1); mrWav1 = -mrWav1 / std(mrWav1(:));
% 
% % spatial convolution
% mrWav2 = zeros(size(mrWav), 'single');
% % miSite = findNearSites(P.mrSiteXY, 1.5);
% % miSite = findNearSites(P.mrSiteXY, 2.5);
% % vrSum = single([4 2 1 1]'/8);
% % vrSum = 
% for iSite = 1:nSites    
% %     mrWav2(:,iSite) = mrWav1(:,miSite(:,iSite)) * vrSum;
% %     mrWav2(:,iSite) = mean(mrWav1(:, P.miSites(:,iSite)), 2);
%     mrWav2(:,iSite) = mean(mrWav(:, P.miSites(:,iSite)), 2);
% end
% mrWav2 = mrWav2 / std(mrWav2(:));
% 
% % plot
% figure; 
% subplot 121; 
% imagesc(mrWav(1:2000,:), [-10 10]);
% hold on; plot(Sevt.viSite, Sevt.viSpk, 'ro');
% title('Filtered');
% 
% subplot 122; 
% imagesc(mrWav2(1:2000,:), [-10 10]);
% hold on; plot(Sevt.viSite, Sevt.viSpk, 'ro');
% title('Filtered-PV1');

% vrWav0 = vrWav0 / std(vrWav0);
% vrWav1 = -vrWav1 / std(vrWav1);
% figure; hold on; 
% plot(vrWav0(1:2000));
% plot(vrWav1(1:2000)); 


function vrWav1 = matched_filt_(vrWav1, vrPv, spkLim)
% vrPv must be filpped already

if isempty(vrPv), return; end
vrWav1 = conv(single(vrWav1), single(vrPv));
% vrWav1 = vrWav1(spkLim(2)+1:end+spkLim(1)); %shift determined by peak loc. shift spike time instead
end %func

% function viSite = nearest_pair_(mrPosXY, mrSiteXY)
% % mrPosXY: nSpk x 2
% % mrR2 = pdist2(mrPosXY, mrSiteXY, 'euclidean'); %use faster algorithm
% mrR2 = eucldist2_(mrSiteXY, mrPosXY);
% [~, viSite] = min(mrR2);
% end %func


function [vrY1, vrY2] = centroid_2col_(mrVpp, mrSiteXY)
% return two column estimate
% [vrX, vrY] = centroid_mr_(mrVpp, vrPos)
% [mrXY] = centroid_mr_(mrVpp, vrPos)
% mrVpp: nSites x nSpk
vl_x0 = mrSiteXY(:,1) == 0;

mrVpp_sq = mrVpp .^ 2;
vrVpp_sq_sum1 = sum(mrVpp_sq(vl_x0,:));
vrVpp_sq_sum2 = sum(mrVpp_sq(~vl_x0,:));
vrY1 = sum(bsxfun(@times, mrVpp_sq(vl_x0,:), mrSiteXY(vl_x0,2))) ./  vrVpp_sq_sum1;
vrY2 = sum(bsxfun(@times, mrVpp_sq(~vl_x0,:), mrSiteXY(~vl_x0,2))) ./  vrVpp_sq_sum2;
if nargout==1
    vrY = [vrY1(:), vrY2(:)];
end
end %func


function [vrY, vrY2] = centroid_mr_(mrVpp, vrYe, mode1)
% [vrX, vrY] = centroid_mr_(mrVpp, vrPos)
% [mrXY] = centroid_mr_(mrVpp, vrPos)
% mrVpp: nSites x nSpk
if nargin<3, mode1 = 1; end
if isrow(vrYe), vrYe=vrYe'; end
% mrVpp_sq = abs(mrVpp);
switch mode1    
    case 1
        mrVpp = abs(mrVpp);
    case 2
        mrVpp = mrVpp.^2; 
    case 3
        mrVpp = sqrt(abs(mrVpp)); 
    case 4
        mrVpp = mrVpp.^2; 
        mrVpp = bsxfun(@minus, mrVpp, min(mrVpp));   
end


vrVpp_sum = sum(mrVpp);
% vrX = sum(bsxfun(@times, mrVpp_sq, mrSiteXY(:,1))) ./  vrVpp_sq_sum;
vrY = sum(bsxfun(@times, mrVpp, vrYe)) ./  vrVpp_sum;
if nargout>=2
    vrY2 = sum(bsxfun(@times, mrVpp, vrYe.^2)) ./  vrVpp_sum;
    vrY2 = sqrt(abs(vrY2 - vrY.^2));
end
% if nargout==1
%     vrX = [vrX(:), vrY(:)];
% end
end %func


function mrR2 = eucldist2_(mrX, mrY)
mrR2 = zeros(size(mrX,1), size(mrY,1));
for i=1:size(mrX,2)
    mrR2 = mrR2 + bsxfun(@minus, mrX(:,i), mrY(:,i)').^2;
end
end


function [cmrFet_site, cviFet_site, miSites_site] = pca_project_site_(mnWav, Sevt, P)
% organized by sites, and each containing spikes centered at sites
fprintf('PCA projection by site (fPcaDetect=1) '); t1=tic;
[nSamples, nSites] = size(mnWav);
% nSpk = numel(Sevt.viSpk);
% nPc = size(P.mrPv,2);
% trFet = zeros([nPc, nSpk, nSites], 'single');
mrPv_tr = Sevt.mrPv';
% mrPv = Sevt.mrPv + eps('single');
% gmrPv_tr = gpuArray(single(P.mrPv'));
% gviSpk = gpuArray(Sevt.viSpk);
[cviFet_site, cmrFet_site] = deal(cell(1, nSites));
nSamples1 = diff(P.spkLim) + 1;

% miSites_site = findNearSites(P.mrSiteXY, P.maxSite+1, P.viSiteZero);
miSites_site = P.miSites;
nSites1 = size(miSites_site,1);
for iSite = 1:nSites
    viSite1 = miSites_site(:, iSite);
    vrPosXY1 = P.mrSiteXY(viSite1,:);
    viSpk1 = find(Sevt.viSite == iSite);
    if isempty(viSpk1), continue; end
    cviFet_site{iSite} = viSpk1;
    viTime1 = Sevt.viSpk(viSpk1);
    mrWav_spk1 = single(mr2tr2(mnWav, viTime1, P.spkLim, 0, viSite1));    
    nSpk1 = numel(viSpk1);
    if P.fMinNorm_wav       
        vrWav_min = std(single(squeeze(mrWav_spk1(:,1,:))));
        mrWav_spk1 = reshape(mrWav_spk1, [], nSpk1);
%         vrWav_min = mrWav_spk1(1-P.spkLim(1), :);
%         viFix = find(vrWav_min==0);
%         vrWav_min(viFix) = min(mrWav_spk1(:,viFix));
%         vrWav_min = min(mrWav_spk1);
%         vrWav_min = single(Sevt.vrSpk(viSpk1));
        mrWav_spk1 = bsxfun(@rdivide, mrWav_spk1, vrWav_min(:)');
    end        
    
        % determine pca across channel
%     viSpk_ap = 1-P.spkLim(1) + (-5:5);
%     mrWav_spk2 = mrWav_spk1(viSpk_ap,:,:);
%     mrWav_spk2 = meanSubt_tr_(mrWav_spk2);
    mrWav_spk1 = meanSubt_tr_(mrWav_spk1, 0);
%     [mrU1, mrC1, vrD1] = trSpk2svd1_(mrWav_spk1, P);  
%     mrPos1 = centroid_mi(mrC1, vrPosXY1);
%     mrPos1 = centroid_vpp(mrWav_spk1, vrPosXY1);
%     [a,b,c]=pca(reshape(mrPos1, size(mrPos1,1), [])');
%     [viShift, mrU2] = align_mr(mrU1, 1-P.spkLim(1) + (-5:5)); %align to mean. subpixel shift after?    
%     cmrFet_site{iSite} = pca(mrU2, 'NumComponents', nSites1)';
%     cmrFet_site{iSite} = mrC1;
%     cmrFet_site{iSite} = reshape((mrPv \ mrWav_spk1)', [], numel(viSpk1));     
    cmrFet_site{iSite} = reshape(mrPv_tr * mrWav_spk1, [], nSpk1);    
    fprintf('.');
end
% trFet = permute(trFet, [1 3 2]);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [mrU1, mrC1, vrD1] = trSpk2svd1_(trSpk1, P)
%     [mrWav_spk2, viShift] = interp_align(mrWav_spk1, P);
%     mrWav_spk1 = reshape(mrWav_spk1, size(mrWav_spk1,1), []); 
%     [U,S,V] = svd(mrWav_spk1 * mrWav_spk1');
%     d=diag(S);
nSpk1 = size(trSpk1,3);
%align at center of the principal     
for iSpk1=1:nSpk1
    mrSpk1 = trSpk1(:,:,iSpk1);
    [U,S,V] = svd(mrSpk1 * mrSpk1', 0);
    if iSpk1==1
        mrU1 = zeros(size(U,1), nSpk1, 'single');
    end
    U1 = U(:,1); %alrady normalized
    U2 = U(:,2); %alrady normalized
    U3 = U(:,3); %alrady normalized
    if U1(1-P.spkLim(1)) > 0, U1 = -U1; end
    if U2(1-P.spkLim(1)) > 0, U2 = -U2; end
    if U3(1-P.spkLim(1)) > 0, U3 = -U3; end
    vrC1 = toVec([U1, U2, U3]' * mrSpk1);
%     vrC1 = U1' * mrSpk1;
    if iSpk1==1
        mrC1 = zeros(numel(vrC1), nSpk1, 'like', trSpk1);
        vrD1 = zeros(nSpk1, 1, 'like', trSpk1);
    end
    mrC1(:,iSpk1) = vrC1;
    mrU1(:,iSpk1) = U1;
%         mrU2(:,iSpk1) = U(:,2);
    d=diag(S);
    vrD1(iSpk1) = d(1)/sum(d);
end
% viFlip1 = find(mrU1(1-P.spkLim(1),:)>0);
%     viFlip1 = find(mrU1(6,:)>0);
% if ~isempty(viFlip1)
%     mrU1(:,viFlip1) = -mrU1(:,viFlip1);   
% end
end
    

function mrWav_spk1 = meanSubt_tr_(mrWav_spk1, fReshape)
if nargin<2, fReshape=1; end
dimm1 = size(mrWav_spk1);
mrWav_spk1 = single(reshape(mrWav_spk1, dimm1(1), []));
mrWav_spk1 = bsxfun(@minus, mrWav_spk1, mean(mrWav_spk1)); %mean subt        
if fReshape
    mrWav_spk1 = reshape(mrWav_spk1, dimm1);
end
end %func


function [trFet, miSites_fet] = matrix_project_(mnWav, Sevt, mrPv, P)
fMeanSpk = 1; %1 for mean subtract, 2 for z score
gmrPv_tr = gpuArray(mrPv)';

[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);
miSites_fet = findNearSites(P.mrSiteXY, P.maxSite*2+1, P.viSiteZero);

nSites1 = size(miSites_fet, 1);
nSamples1 = diff(P.spkLim) + 1;
trFet = zeros([P.nPcPerChan, nSites1, nSpk], 'single');
for iSite = 1:nSites
    viiSpk1 = find(Sevt.viSite==iSite);
    nSpk1 = numel(viiSpk1);
    if nSpk1==0, continue; end
    viSites1 = miSites_fet(:,iSite);
    gtrFet1 = gpuArray(mr2tr3(mnWav, P.spkLim, Sevt.viSpk(viiSpk1), viSites1));    
    if P.fMeanSite == 1 %mean across site subtraction
        gtrFet1 = trSpk_mean_site(gtrFet1, 3);
    elseif P.fMeanSite == 2
        gtrFet1 = trSpk_med_site(gtrFet1, 3);
    end
    gtrFet1 = single(reshape(gtrFet1, nSamples1, []));
    if fMeanSpk == 1
        gtrFet1 = bsxfun(@minus, gtrFet1, mean(gtrFet1)); %mean subtract
    elseif fMeanSpk == 2
        gtrFet1 = zscore(gtrFet1);
    end
    gtrFet1 = gmrPv_tr * gtrFet1;   % project 
    gtrFet1 = reshape(gtrFet1, [], nSpk1, nSites1);
    trFet(:,:,viiSpk1) = gather(permute(gtrFet1, [1, 3, 2]));
    fprintf('.');
end
% trFet = permute(trFet, [1 3 2]);
end %func


function [trFet, miSites_fet] = diff_project_(mnWav, Sevt, P) % P.mrPv used
[trFet, miSites_fet] = diff_project(mnWav, Sevt, P);
end %func


function [trFet] = fet_tetrode_(mnWav, Sevt, P)
fprintf('Tetrode feature \n'); t1=tic;
[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);


mrSpk = single(permute(mr2tr3(mnWav, [-12 12], Sevt.viSpk), [1,3,2]));
dimm1 = size(mrSpk);
% mrVpp = squeeze(max(mrSpk)-min(mrSpk));
mrSpk = reshape(mrSpk, size(mrSpk,1), []);
[trFet,~,~]=pca(mrSpk, 'NumComponents', P.nPcPerChan);
trFet = permute(reshape(trFet, [dimm1(2:3), P.nPcPerChan]), [3,1,2]);
% trFet = reshape(trFet, [], dimm1(3));
% gtrWav1 = meanSubt(gtrWav1);
% dimm1 = size(gtrWav1);
% gtrWav1 = zscore_tr(gtrWav1); 
% 
% gtrWav1.*gtrWav1
% 
% fprintf('\n\ttook %0.1fs\n', toc(t1));
% trFet = gather(trFet);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [trFet, miSites_fet] = fet_spacetime_(mnWav, Sevt, P) % P.mrPv used
% output features that are stationary. detrend all the non-stationary
% signal.
% detrend using y depth first and apply trending after sorting by depth
% determine x,y,z,t and properly normalize by max abs dist
% viRange_fet = [];
% if isfield(P, 'spkLim_ms_fet')
%     if ~isempty(P.spkLim_ms_fet)
%         spkLim_fet = round(P.spkLim_ms_fet * P.sRateHz / 1000);
%         viRange_fet = spkLim_fet + 1 - P.spkLim(1);
%     end
% end
% fMaxNorm = 1;

fprintf('Spacetime feature \n'); t1=tic;
[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);
nSites1 = min(P.maxSite_fet*2+1, nSites); % 2 + 4*n
fPca = 0;
nRef = 0;
miSites_fet = findNearSites(P.mrSiteXY, P.maxSite_fet + nRef/2, P.viSiteZero);
% [vrX_cm, vrY_cm, vrA_cm] = deal(zeros(nSpk, 1, 'single', 'gpuArray'));
trFet = [];
% nRms_filt = round(P.rms_filt_ms*P.sRateHz/1000);
% spkLim_fet = -nRms_filt:nRms_filt;
% spkLim_fet = [-P.nDiff_filt,P.nDiff_filt*2];  spkLim_fet2 = [P.nDiff_filt,P.nDiff_filt*3];
% spkLim_fet = [-P.nDiff_filt,0];  spkLim_fet2 = [0,P.nDiff_filt];

% spkLim_fet = [-12,12];
% spkLim_fet = round(P.spkLim_fet_ms * P.sRateHz / 1000);
%vrWin = hann(numel(spkLim_fet));
% vrWin = [];
% viRange1 = 1-P.spkLim(1) + (-nRms_filt:nRms_filt);
% viFft = 2:round(nRms_filt/2);
% abssum = @(x)sum(abs(x));
% mr_fft = dftmtx(numel(spkLim_fet));
if fPca
    viSubs = subsample_vr(1:numel(Sevt.viSpk), 10000);
    tr2Vpp(mr2tr4(mnWav, P.spkLim_fet, Sevt.viSpk(viSubs), Sevt.viSite(viSubs)));
end
for iSite = 1:nSites
    try
        viiSpk1 = find(Sevt.viSite==iSite);    
        nSpk1 = numel(viiSpk1);
        if nSpk1==0, continue; end
        viSites0 = miSites_fet(:, iSite);
        viSites1 = viSites0(1:nSites1);
        gmrXYe = (single(P.mrSiteXY(viSites1,:)));
        viTime1 = Sevt.viSpk(viiSpk1);
        gtrWav1 = mr2tr3(mnWav, P.spkLim_fet, viTime1, viSites1);
        gtrWav1 = single(gtrWav1);
        if nRef>0
            gmrWav1_ref = mean(gtrWav1(:,:,end-nRef+1:end), 3); %surrounding ground
            gtrWav1 = gtrWav1(:,:,1:end-nRef) - repmat(gmrWav1_ref, [1,1,nSites1]); %csd
        end                
        vlLeft = gmrXYe(:,1) < mean(gmrXYe(:,1));
        [mrVpp1, mrVpp2, mrVpp3] = tr2Vpp(gtrWav1, P);
        mrFet1 = permute(cat(3, mrVpp1, mrVpp2, mrVpp3), [1,3,2]);
        vyPos_spk1 = centroid_mr_(mrVpp1, gmrXYe(:,2), 2);
        vxPos_spk1 = centroid_mr_(mrVpp1, gmrXYe(:,1), 2);
        if 1
            vrSiteOrder = gmrXYe(:,2) + gmrXYe(:,1)*max(gmrXYe(:,2));
            [~, viSiteOrder] = sort(vrSiteOrder);
            mrFet1 = mrFet1(viSiteOrder,:,:);
        elseif 0
            [~, viSiteOrder] = sort(viSites1, 'ascend');
            mrFet1(viSiteOrder,:,:) = mrFet1; %reshape(mrFet1(viSiteOrder,:,:), [], nSpk1);
        else
            mrSiteDist_spk1 = bsxfun(@minus, gmrXYe(:,2), vyPos_spk1).^2 + bsxfun(@minus, gmrXYe(:,1), vxPos_spk1).^2;
            [mrSiteDist_spk1, miSrt1] = sort(mrSiteDist_spk1, 1, 'ascend');
            for iSpk1=1:nSpk1
                mrFet1(:,:,iSpk1) = mrFet1(miSrt1(:,iSpk1),:,iSpk1);
            end            
        end
        mrFet1 = reshape(mrFet1, [], nSpk1);
%         mrFet1 = signsqrt(mrFet1);            
%         mrFet1 = (abs(mrFet1)); %square weighing
        mrFet1 = [vyPos_spk1; vxPos_spk1; mrFet1]; %incl position 
    %      mrFet1 = [centroid_mr_(mrVpp1, gmrXYe(:,2), 2) + centroid_mr_(mrVpp1, gmrXYe(:,1), 2); mrFet1];
         if isempty(trFet)
            trFet = zeros([size(mrFet1,1), nSpk], 'single', 'gpuArray'); 
         end
         trFet(:, viiSpk1) = mrFet1;
        fprintf('.');
    catch
        disp(lasterr());
    end
end
fprintf('\n\ttook %0.1fs\n', toc(t1));
trFet = gather(trFet);

%----------------------------------
% drift correction
viSpk = single(Sevt.viSpk);
vrY_cm = trFet(1,:)';
vrA_cm = mean(abs(trFet(3:end,:)),1); %exp(trFet(2,:))' + exp(trFet(end,:))';

if ~isempty(P.nT_drift)
    nT_cm = P.nT_drift;
elseif ~isempty(P.tbin_drift)
    t_dur = (max(viSpk)-min(viSpk)) / P.sRateHz;
    nT_cm = floor(t_dur / P.tbin_drift);
else
    nT_cm = []; %do not perform drift correction
end
if ~isempty(nT_cm)
    nY_cm = (max(P.mrSiteXY(:,2)) - min(P.mrSiteXY(:,2))) / P.ybin_drift; %spatial bin
end
try
    if ~isempty(nT_cm)
        fprintf('Drift correcting\n'); t_drift=tic;
        fPlot = 1;
        while 1            
            vrY_cm_pre = vrY_cm;
            [vrY_cm, ~, ~] = unshift_spikes(viSpk, vrY_cm, vrA_cm, nT_cm, nY_cm, 'nneigh', 2, 'fPlot', fPlot);
            fPlot = 0;
            nT_cm = round(nT_cm/2);
            if nT_cm < 2 
                unshift_spikes(viSpk, vrY_cm_pre, vrA_cm, nT_cm*2, nY_cm, 'nneigh', 2, 'fPlot', 1);
                break; 
            end                        
            fprintf('.');            
        end
        fprintf('\n\ttook %0.1fs\n', toc(t_drift));
        trFet(1,:) = vrY_cm; %update the cm
    end
catch
    disp(lasterr());
end
trFet = reshape(trFet, [1, size(trFet)]);
miSites_fet = []; %common featureset for all sites
end %func


function [trFet, miSites_fet] = fet_spacetime_old_(mnWav, Sevt, P) % P.mrPv used


fprintf('Spacetime feature\n'); t1=tic;
[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);
% trFet = zeros([(1 + numel(viDelay)*2), nSpk], 'single', 'gpuArray');
% mrFet(:,end) = single(Sevt.viSpk); %time index
nSites1 = P.maxSite*2+1; % 2 + 4*n
P.fMeanSite = 0;
nRef = 0;
% miSites_fet = findNearSites(P.mrSiteXY, nSites1*2+1, P.viSiteZero);
miSites_fet = findNearSites(P.mrSiteXY, P.maxSite + nRef/2, P.viSiteZero);
% [mrA_site] = zeros(nSites, 2, 'single', 'gpuArray');
% figure; hold on;
% [cmrFet_sub_site, cviSpk_sub_site, cvrDk_sub_site] = deal(cell(1, nSites));
% , 
[vrX_cm, vrY_cm, vrA_cm] = deal(zeros(nSpk, 1, 'single', 'gpuArray'));
% viNneigh = zeros(nSpk, 1, 'int32');
% determine amplitudes and subsample by channel
trFet = [];
for iSite = 1:nSites
    try
    viiSpk1 = find(Sevt.viSite==iSite);    
    nSpk1 = numel(viiSpk1);
    if nSpk1==0, continue; end
    viSites0 = miSites_fet(:, iSite);
    viSites1 = viSites0(1:nSites1);
    gmrXYe = gpuArray(single(P.mrSiteXY(viSites1,:)));
%     [gmrV_diff, gmrV] = autoproj_(mnWav, Sevt.viSpk(viiSpk1), viSites1, P);
    viTime1 = Sevt.viSpk(viiSpk1);
%     gtrWav1 = mr2tr3(mnWav, P.spkLim, viTime1, viSites1);
%     gmrV = tr2pc1(gtrWav1);
%     gtrWav1 = cumsum(gtrWav1, 1);
%     gmrV = shiftdim(std(single(gtrWav1),1,1), 1);
%     [gmrV] = autoproj_(mnWav, viTime1, viSites1, viDelay, P);

    %[mrMin1, mrMax1] = peak_maxsite(mnWav, viTime1, viSites1, P); %returns a cell
    gtrWav1 = mr2tr3(mnWav, P.spkLim, viTime1, viSites0);
    gtrWav1 = single(gtrWav1);
%     gtrWav1 = gtrWav1 - repmat(mean(gtrWav1,1), [size(gtrWav1,1),1,1]);
    
    % subtract surrounding
%     gtrWav1 = single(gtrWav1);
        
    if nRef>0
        gmrWav1_ref = mean(gtrWav1(:,:,end-nRef+1:end), 3); %surrounding ground
        gtrWav1 = gtrWav1(:,:,1:end-nRef) - repmat(gmrWav1_ref, [1,1,nSites1]); %csd
    end
%     gtrWav1 = cumsum(gtrWav1,1);
%     gmrXYe(1,:) = [];
    
%     mrVpp = single(squeeze(max(gtrWav1)-min(gtrWav1))'); %peak to peak works better than std
%     mrVpp = single(squeeze(min(gtrWav1))'); %peak to peak works better than std
    it0 = 1 - P.spkLim(1);
    vlXa_left = (gmrXYe(:,1) <= mean(gmrXYe(:,1),1));
    
    mrVmax = squeeze(max(gtrWav1(it0+1:end,:,:)))';
    mrVmin = squeeze(min(gtrWav1(1:it0,:,:)))';
    mrVpp = single(abs(mrVmax-mrVmin));


    vrSiteOrder = gmrXYe(:,2) + gmrXYe(:,1)*10;
    [~, viSiteOrder] = sort(vrSiteOrder);
    mrFet1 = mrVpp(viSiteOrder, :);
    
    mrFet1 = [centroid_mr_(mrVpp, gmrXYe(:,2), 2); sqrt(abs(mrFet1))];
    
%     mrFet1 = [centroid_mr_(mrVpp', gmrXYe(:,2), 2); ...
%               sqrt(elec_quadrant(mrVpp, gmrXYe))];
%     mrFet1 = [centroid_mr_(mrVpp', gmrXYe(:,2), 2); ...
%               sqrt(abs(mrVpp'))];
%     mrFet1 = [centroid_mr_(gmrV{1}', gmrXYe(:,2), 2); ...
%                          centroid_mr_(gmrV{1}', gmrXYe(:,1), 2); ...
%                          std(gmrV{1},1,2)'; std(gmrV{2},1,2)'];
     if isempty(trFet)
        trFet = zeros([size(mrFet1,1), nSpk], 'single', 'gpuArray'); 
     end
     trFet(:, viiSpk1) = mrFet1;
%                          elec_quadrant(gmrV{1}, gmrXYe)];
    fprintf('.');
    catch
        disp(lasterr());
    end

end
fprintf('\n\ttook %0.1fs\n', toc(t1));
trFet = gather(trFet);

%----------------------------------
% drift correction
viSpk = single(Sevt.viSpk);
vrY_cm = trFet(1,:)';
vrA_cm = mean(abs(trFet(2:end,:)),1); %exp(trFet(2,:))' + exp(trFet(end,:))';

if ~isempty(P.nT_drift)
    nT_cm = P.nT_drift;
elseif ~isempty(P.tbin_drift)
    t_dur = (max(viSpk)-min(viSpk)) / P.sRateHz;
    nT_cm = round(t_dur / P.tbin_drift);
else
    nT_cm = []; %do not perform drift correction
end
nY_cm = (max(P.mrSiteXY(:,2)) - min(P.mrSiteXY(:,2))) / P.ybin_drift; %spatial bin

try
    disp('drift correcting');
    [vrY_cm_corr, ~, ~] = unshift_spikes(viSpk, vrY_cm, vrA_cm, nT_cm, nY_cm, 'nneigh', 2);
    [vrY_cm_corr, ~, ~] = unshift_spikes(viSpk, vrY_cm_corr, vrA_cm, round(nT_cm/2), nY_cm, 'nneigh', 2);
    [vrY_cm_corr, ~, ~] = unshift_spikes(viSpk, vrY_cm_corr, vrA_cm, round(nT_cm/4), nY_cm, 'nneigh', 2);
    [vrY_cm_corr, ~, ~] = unshift_spikes(viSpk, vrY_cm_corr, vrA_cm, round(nT_cm/8), nY_cm, 'nneigh', 2);
    trFet(1,:) = vrY_cm_corr; %update the cm
catch
    disp(lasterr());
end
% [vrY_cm_corr, vrA_cm] = multifun(@gather, vrY_cm_corr, vrA_cm);

% mrFet = [vrY_cm_corr(:), mrFet];
trFet = reshape(trFet, [1, size(trFet)]);
% mrFet(:,2) = correct_depth(vrT_sub, vrY_cm(viSpk_sub), vrDepthY, P_shift);
miSites_fet = []; %common featureset for all sites
end %func


function gmrFet1 = autoproj_(mnWav, viTime1, viSites1, viDelay, P)
% 
% nSites1 = size(miSites_fet, 1);
% nSamples1 = diff(P.spkLim) + 1;

% viDelay = [-2,0,2];
fSiteAve = 3; %0: use largest site as a template
fIntegrate = 1;
fDiffMax = 0;

% determine spike waveform
gtrWav1 = mr2tr3(mnWav, P.spkLim, viTime1, viSites1);  
gtrWav1 = single(gpuArray(gtrWav1));
% if P.fMeanSite == 1 %mean across site subtraction
%     gtrWav1 = trSpk_mean_site(gtrWav1, 3);
% elseif P.fMeanSite == 2
%     gtrWav1 = trSpk_med_site(gtrWav1, 3);
% end
if fIntegrate
%     gtrWav1 = meanSubt(gtrWav1);
    gtrWav1 = cumsum(gtrWav1, 1);
end
% deterime sites toproject
switch fSiteAve
    case 1
        for iSiteAve = 1:3
            if iSiteAve==1
                gmrWav1 = zscore(gtrWav1(:,:,iSiteAve));
            else
                gmrWav1 = gmrWav1 + zscore(gtrWav1(:,:,iSiteAve));
            end
        end
        gmrWav1 = gmrWav1 / 3;
    case 2
        gmrWav1 = mean(gtrWav1, 3);
    case 3
        gtrWav2 = repmat(gtrWav1(:,:,1), [1,1,size(gtrWav1,3)-1]) - gtrWav1(:,:,2:end);
        gmrWav1 = mean(gtrWav2, 3);
    case 0
    gmrWav1 = gtrWav1(:,:,1); %first site
end
if fDiffMax
    gtrWav1 = repmat(gmrWav1, [1,1,size(gtrWav1,3)-1]) - gtrWav1(:,:,2:end);
end
% gtrWav1 = permute(gtrWav1, [1,3,2]);        
[gmrFet1, gtrWav1] = project_tr2mr(gtrWav1, gmrWav1, viDelay);
% gmrFet1 = gmrFet1';
% else
%     gmrFet2 = gmrFet2 / sqrt(2); %scaling
% end
% end
end %func


function [trFet, miSites_fet] = autocor_project_(mnWav, Sevt, P) % P.mrPv used
% todo: time delay
fprintf('self corr feature (fPcaDetect=1) '); t1=tic;
fSiteAve = 0;
viShift = [0, 2, -2, 1, -1];
viShift = viShift(1:P.nPcPerChan);

[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);
% miSites_fet = findNearSites(P.mrSiteXY, P.maxSite*2+1, P.viSiteZero);
nSites2 = 2+4*3; % 2 + 4*n
nSites1 = 6;

P.fMeanSite = 0;
miSites_fet2 = findNearSites(P.mrSiteXY, (nSites2-1)/2, P.viSiteZero);
miSites_fet = miSites_fet2(1:nSites1, :);
% nSites1 = size(miSites_fet, 1);
nSamples1 = diff(P.spkLim) + 1;
% trFet = zeros([P.nPcPerChan, nSites1, nSpk], 'single');
trFet = zeros([1, nSites1, nSpk], 'single');
for iSite = 1:nSites
    viiSpk1 = find(Sevt.viSite==iSite);
    nSpk1 = numel(viiSpk1);
    if nSpk1==0, continue; end
    viSites1 = miSites_fet2(:,iSite);
    gtrWav1 = mr2tr3(mnWav, P.spkLim, Sevt.viSpk(viiSpk1), viSites1);  
    gtrWav1 = single(gpuArray(gtrWav1));
    if P.fMeanSite == 1 %mean across site subtraction
        gtrWav1 = trSpk_mean_site(gtrWav1, 3);
    elseif P.fMeanSite == 2
        gtrWav1 = trSpk_med_site(gtrWav1, 3);
    end
%     gmrSpk1_tr = shift_mr(gtrWav1(:,:,1), viShift); %first site    
    if fSiteAve
        for iSiteAve = 1:3
            if iSiteAve==1
                gmrWav1 = zscore(gtrWav1(:,:,iSiteAve));
            else
                gmrWav1 = gmrWav1 + zscore(gtrWav1(:,:,iSiteAve));
            end
        end
        gmrWav1 = gmrWav1 / 3;
    else
        gmrWav1 = gtrWav1(:,:,1);
    end
%     gmrWav11_tr = zscore(gmrWav1') / nSamples1;    

%     gtrWav1 = tr_mean_subt(gtrWav1, 1); % mean subtract    
    gtrWav1 = permute(gtrWav1, [1,3,2]);    
%     if P.nPcPerChan == 1
%         gtrFet1 = zeros([P.nPcPerChan, nSites1, nSpk1], 'single', 'gpuArray');
%         for iSpk1=1:nSpk1 %make it faster.
%             gtrFet1(:,:,iSpk1) = gmrWav11_tr(iSpk1,:) * gtrWav1(:,:,iSpk1);
%         end
%     else
    gmrXYe = P.mrSiteXY(viSites1, :);
        gtrFet1 = zeros([nSites1, nSpk1, P.nPcPerChan], 'single', 'gpuArray');
        for iFet1 = 1:P.nPcPerChan
            gmrWav11 = zscore(shift_mr(gmrWav1, viShift(iFet1))) / nSamples1;
            gmrV1 = proj_spk(gmrWav11, gtrWav1);
            [gmrV1, gmrXYe_a, gmrXYe_b] = ring_ref(gmrV1, gmrXYe, nSites1);
            gtrFet1(:,:,iFet1) = gmrV1;
        end        
        gtrFet1 = permute(gtrFet1, [3,1,2]);
%         for iSpk1=1:nSpk1 %make it faster.
%             gtrFet1(:,:,iSpk1) = zscore(shift_mr(gmrWav1(:,iSpk1), viShift))' * ...
%                 gtrWav1(:,:,iSpk1);
%         end
%     end
    trFet(:,:,viiSpk1) = gather(gtrFet1);
    fprintf('.');
end

fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [trFet, miSites_fet] = moment_project_(mnWav, Sevt, P) % P.mrPv used
fprintf('moment feature (fPcaDetect=1) '); t1=tic;
vrPow = [1, 2, 4, .5, .25]; %powers
vrPow = vrPow(1:P.nPcPerChan);

% build a time-shifted matrix
vrPv_abs = abs(Sevt.vrPv_mean(:));
vrPv_sign = sign(Sevt.vrPv_mean(:));
mrPv = zeros(numel(vrPv_abs), numel(vrPow), 'single');
for iPow = 1:numel(vrPow)
    mrPv(:,iPow) = vrPv_sign .* (vrPv_abs .^ vrPow(iPow));
end
mrPv = zscore(mrPv) * P.uV_per_bit / numel(Sevt.vrPv_mean);

[trFet, miSites_fet] = matrix_project_(mnWav, Sevt, mrPv, P);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [trFet, miSites_fet] = xcor_project_(mnWav, Sevt, P) % P.mrPv used
fprintf('xcor feature (fPcaDetect=1) '); t1=tic;
% fMeanSpk = 1; %1 for mean subtract, 2 for z score
fProjMean_shift = [0 1 -1 2 -2];
fProjMean_shift = fProjMean_shift(1:P.nPcPerChan);

% build a time-shifted matrix
vrPv = Sevt.vrPv_mean(:);
mrPv = repmat(vrPv, [1, numel(fProjMean_shift)]);
mrPv = shift_mr(mrPv, fProjMean_shift);
mrPv = zscore(mrPv) * P.uV_per_bit / size(mrPv,1);

[trFet, miSites_fet] = matrix_project_(mnWav, Sevt, mrPv, P);

fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function [trFet, miSites_fet] = pca_project_(mnWav, Sevt, P) % P.mrPv used
fprintf('PCA projection (fPcaDetect=1) '); t1=tic;
[trFet, miSites_fet] = matrix_project_(mnWav, Sevt, Sevt.mrPv, P);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


function trFet = pca_project_all_(mnWav, Sevt, P) % P.mrPv used
fprintf('PCA projection (fPcaDetect=1) '); t1=tic;
[nSamples, nSites] = size(mnWav);
nSpk = numel(Sevt.viSpk);
nPc = size(Sevt.mrPv,2);
trFet = zeros([nPc, nSpk, nSites], 'single');

gmrPv_tr = gpuArray(Sevt.mrPv');
% gmrPv = gpuArray(Sevt.mrPv + eps('single'));
gviSpk = gpuArray(Sevt.viSpk);
for iSite = 1:nSites
    % project waveforms
%     gvnWav1 = gpuArray(mnWav(:,iSite));
    gmrWav_spk1 = single(vr2mr2(gpuArray(mnWav(:,iSite)), gviSpk, P.spkLim));
%     gmrWav_spk1 = squeeze(mr2tr2(mnWav, Sevt.viSpk, P.spkLim, 0, iSite));
%     gmrWav_spk1 = single(gpuArray(gmrWav_spk1));
    trFet(:,:,iSite) = gather(gmrPv_tr * gmrWav_spk1); 
%     trFet(:,:,iSite) = gather(gmrPv \ gmrWav_spk1)'; 
    fprintf('.');
end
trFet = permute(trFet, [1 3 2]);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


% function Sclu = recenter_Sevt_(Sevt, P)
% % determine spike centroid and shift the center
% % need cmrVpp_site before hand
% nSites = numel(Sevt.cviSpk_site);
% miSite = findNearest(P.mrSiteXY, P.maxSite);
% [cvrPosX_site, cvrPosY_site] = deal(cell(1, nSites));
% 
% for iSite = 1:nSites
%     % Turn centroid into nearest site number
%     mrSiteXY1 = P.mrSiteXY(miSite(:, iSite), :);
%     [cvrPosX_site{iSite}, cvrPosY_site{iSite}] = ...
%         centroid_mr_(Sevt.cmrVpp_site{iSite}, mrSiteXY1);
% end
% Sclu.cvrPosX_site = cvrPosX_site;
% Sclu.cvrPosY_site = cvrPosY_site;
% end %func


function vr = cell2mat_(cvr)
% remove empty
vi = find(cellfun(@(x)~isempty(x), cvr));
vr = cell2mat(cvr(vi));
end %func


function export_workspace_()
% Exports fields S0 to workspace

S0 = get(0, 'UserData');
if isempty(S0), disp('Nothing to export.'); return; end
csNames = fieldnames(S0);
for i=1:numel(csNames)
    if strcmpi(csNames{i}, 'mrWav'), continue; end
    eval(sprintf('%s = S0.%s;', csNames{i}, csNames{i}));
    eval(sprintf('assignWorkspace(%s);', csNames{i}));
end
vcFile_prm = '';
try
    if isfield(S0, 'P'), vcFile_prm = P.vcFile_prm; end
    if isfield(S0, 'Sevt'), vcFile_prm = Sevt.P.vcFile_prm; end
    if isfield(S0, 'Sclu'), vcFile_prm = Sclu.P.vcFile_prm; end
catch
    ;
end
fprintf('Exported %s\n', vcFile_prm);
end %func


function import_clusters_(vcFile_prm)
Sgt = struct();
if ~isempty(vcFile_prm)
    P0 = loadParam_(vcFile_prm);
    [vcDir, ~, ~] = fileparts(P0.vcFile);
else
    vcDir = '';
end
vcFile_time = uiFileDialog_([vcDir, filesep(), '*'], 'PROVIDE SPIKE TIMES');
if isempty(vcFile_time), return; end
if matchFileExt(vcFile_time, '.csv') %catalin in-silico
   % catalin's csv format
    mrCsv = importdata(vcFile_time);
    Sgt.viTime = round(mrCsv(:,1) * 20000);
    Sgt.viClu = mrCsv(:,2);
    if size(mrCsv,2)>=3, Sgt.viSite = mrCsv(:,3); end
elseif matchFileExt(vcFile_time, '.bin')
    % Adam Kampff intracell
    [vcDir, vcFile_bin, ~] = fileparts(vcFile_time);
    P = struct('nChans', 8, 'uV_per_bit', 1.5259, 'thresh', [], 'iChan', 1, 'vcDataType', 'uint16', 'sRateHz', 30000, 'fIntra', 1, 'nSpk', []);
    if issubstr(vcDir, '2015_09_03_Pair_9_0'), P.iChan = 1; P.vcDataType = 'uint16'; P.fIntra=0;  P.thresh = 150;
    elseif issubstr(vcDir, '2014_11_25_Pair_3_0'), P.iChan = 1; %P.nSpk = 348;
    elseif issubstr(vcDir, '2014_03_26_Pair_2_0'), P.iChan=2; P.vcDataType='int32';
    elseif issubstr(vcDir, '2014_03_26_Pair_2_1'), P.iChan=2; P.vcDataType='int32';
    end    
    %fid = fopen(vcFile_time); mrWav = fread(fid, inf, sprintf('%s=>single',P.vcDataType)); fclose(fid);
    mrWav = loadWavFile(vcFile_time, P.nChans, inf, P.vcDataType);
    vrWav = single(mrWav(:,P.iChan)) * P.uV_per_bit;
    if P.fIntra
        %vrWav = diff_filt(vrWav); 
    else        
        vrWav = filtfilt_chain(vrWav, 'sRateHz', P.sRateHz, 'freqLim', [300 3000]);        
        vrWav = fft_clean(vrWav, 10);
        vrWav = -(vrWav); 
    end
   % vrWav = filtfilt_chain(vrWav, 'freqLim', [300 3000], 'filtOrder', 3, 'sRateHz', P.sRateHz);    
    figure; hold on; plot(vrWav);
    if isempty(P.thresh)
        P.thresh = (max(vrWav) + median(vrWav))/2;
    end
    [viTime, vrSpk] = findPeak(vrWav, P.thresh);
    [viTime, vrSpk] = spike_refrac(viTime, vrSpk, P0.spkRefrac); %remove refractory spikes
    if ~isempty(P.nSpk) && numel(viTime) > P.nSpk
        [~, viSrt] = sort(vrSpk, 'descend'); 
        viTime = viTime(viSrt(1:P.nSpk)); 
    end
    plot(viTime, vrWav(viTime), 'ro'); title(sprintf('# spikes: %d', numel(viTime)));
    plot(viTime([1,end]), [1,1]*P.thresh, 'r-');
    set(gcf,'Name', vcFile_time, 'color','w'); grid on;
    
    Sgt.viClu = ones(size(viTime));
    Sgt.viTime = viTime;
elseif matchFileExt(vcFile_time, '.mat')
    % catalin raster
    S = load(vcFile_time);    
    vnSpk = cellfun(@numel, S.a);    
    
    vcFile_csv = uiFileDialog_([vcDir, filesep(), '*'], 'Provide cell index file (cancel to skip)');
    if ~isempty(vcFile_csv)
        viCell = importdata(vcFile_csv);
    	S.a = S.a(viCell);
        vnSpk = vnSpk(viCell);
    end
    Sgt.viClu = cell2mat(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0));
    Sgt.viTime = cell2mat(S.a) * 20;
else
    % nick steinmetz txt format
    Sgt.viTime = importdata(vcFile_time);

    vcFile_clu = uiFileDialog_([vcDir, filesep(), '*'], 'PROVIDE CLUSTER NUMBERS');
    if isempty(vcFile_clu), return; end
    Sgt.viClu = importdata(vcFile_clu);

    csAns = inputdlg('Cluster number to track (provide a cell array, leave blank to skip)', '');
    if ~isempty(csAns{1})
        eval(sprintf('cviClu_gt = %s;', csAns{1}));
        [cviTime, cviClu] = deal(cell(size(cviClu_gt')));
        for iGt=1:numel(cviClu_gt)
            viClu_gt1 = cviClu_gt{iGt};
    %         vi_gt1 = find(ismember(Sgt.viClu, viClu_gt1));
            cviTime{iGt} = Sgt.viTime(ismember(Sgt.viClu, viClu_gt1));
            cviClu{iGt} = iGt * ones(size(cviTime{iGt}));
        end
        Sgt.viClu = cell2mat(cviClu);
        Sgt.viTime = cell2mat(cviTime);
    end
    if numel(Sgt.viClu) > numel(Sgt.viTime)
        Sgt.viClu(1) = []; 
    end
end

% determine triggered using unfiltered traces
% try
%     if matchFileExt(vcFile_prm, '.prm')
%         P = loadParam_(vcFile_prm);    
%         P.viChan = P.viSite2Chan;
%         P.freqLim = [];
%         P.nDiff_filt = 0;
%         P.nSamples_load = 1e7; %complications with transpose if divided
%         mrWav = import_whisper(P.vcFile, P); %raw channel read
%         Sgt.tmrCluWav_mu = Sclu_meanWav(Sgt, mrWav, P);
%         Sgt.mrVpp_clu = squeeze(max(Sgt.tmrCluWav_mu)-min(Sgt.tmrCluWav_mu));
%         [Sgt.vrVpp_clu, Sgt.viSite_clu] = max(Sgt.mrVpp_clu, [], 1);
%     end
% catch
%     disp(lasterr);
% end
% try
    vcFile_gt = subsFileExt(vcFile_prm, '_gt.mat');
% % catch
%     vcFile_gt = subsFileExt(vcFile_prm, '_gt.mat');
% end
vcAns = inputdlg('Ground-truth file name', 'Confirmation', 1, {vcFile_gt});
if isempty(vcAns), return; end
vcFile_gt = vcAns{1};
save(vcFile_gt, '-struct', 'Sgt');
% msgbox(['Exported to ', vcFile_gt]);
disp(['Exported to ', vcFile_gt]);
end %funce


function import_openephys_(vcFile)
% import openephys format
[vcDir,~,~]=fileparts(vcFile);
csFiles = dir([vcDir, filesep(), '*.continuous']);
if isempty(csFiles)
    fprintf(2, 'No files found.\n'); 
    return; 
end
nChans = numel(csFiles);

% create a .bin file (fTranspose = 0)
fprintf('Converting to WHISPER format (.bin and .meta)\n\t'); t1=tic;
mnWav = [];
for iChan=1:nChans
    vcFile1 = sprintf('%s\\100_CH%d.continuous', vcDir, iChan);
    if ~exist(vcFile1, 'file'), continue; end
    if isempty(mnWav)
        [vnWav1, ~, Smeta] = load_open_ephys_data_faster(vcFile1, 'unscaledInt16');
        mnWav = zeros(numel(vnWav1), nChans, 'like', vnWav1);
        mnWav(:,iChan) = vnWav1;
    else
        mnWav(:,iChan) = load_open_ephys_data_faster(vcFile1, 'unscaledInt16');
    end
    fprintf('.');
end
fprintf('\n\ttook %0.1fs.\n', toc(t1));

% output bin
mnWav=mnWav';
vcFile_bin = [vcDir, filesep(), 'openephys.bin'];
fid = fopen(vcFile_bin, 'W');
fwrite(fid, mnWav, 'int16');
fclose(fid);

% write .meta file
vcFile_meta = [vcDir, filesep(), 'openephys.meta'];
fid = fopen(vcFile_meta, 'W');
fprintf(fid, 'niMNGain=200\n');
fprintf(fid, 'niSampRate=30000\n'); %intan hardware default. in Smeta.header
fprintf(fid, 'niAiRangeMax=1.278\n'); %intan hardware default. in Smeta.header
fprintf(fid, 'niAiRangeMin=-1.278\n'); %intan hardware default. in Smeta.header
fprintf(fid, 'nSavedChans=%d\n', nChans); %intan hardware default. in Smeta.header
fclose(fid);
% write meta file and bin file

fprintf('Exported to %s, %s\n', vcFile_bin, vcFile_meta);
end %func


function import_tsf_(vcFile)
% import tsf format (test spike file)
% create a .bin file (fTranspose = 0)
fprintf('Converting to WHISPER format (.bin and .meta)\n\t'); t1=tic;
[mnWav, Sfile] = importTSF(vcFile);
nChans = size(mnWav,2);
fprintf('\n\ttook %0.1fs.\n', toc(t1));

% output bin
vcFile_bin = subsFileExt(vcFile, '.bin');
fid = fopen(vcFile_bin, 'W');
fwrite(fid, mnWav', 'int16'); %must transpose
fclose(fid);

% write .meta file
vcFile_meta = subsFileExt(vcFile, '.meta');
fid = fopen(vcFile_meta, 'W');
fprintf(fid, 'niMNGain=200\n');
fprintf(fid, 'niSampRate=%d\n', Sfile.sRateHz); %intan hardware default. in Smeta.header
fprintf(fid, 'niAiRangeMax=0.6554\n'); %intan hardware default. in Smeta.header
fprintf(fid, 'niAiRangeMin=-0.6554\n'); %intan hardware default. in Smeta.header
fprintf(fid, 'nSavedChans=%d\n', Sfile.nChans); %intan hardware default. in Smeta.header
fprintf(fid, 'fileTimeSecs=%f\n', Sfile.n_vd_samples/Sfile.sRateHz); %intan hardware default. in Smeta.header
fprintf(fid, 'fileSizeBytes=%d\n', round(Sfile.n_vd_samples*2*Sfile.nChans)); %intan hardware default. in Smeta.header
fclose(fid);

% write meta file and bin file
assignWorkspace(Sfile);
fprintf('Exported to %s, %s\n', vcFile_bin, vcFile_meta);
end %func


function verify_detection_(vcFile)
% load timing and see how many ground truth clusters are well detected

vcFile_gt = subsFileExt(vcFile, '_gt.mat'); %ground truth file    
Sgt = load_gt_(vcFile_gt);    
if isempty(Sgt), fprintf(2, 'Groundtruth does not exist. Run "jrclust import" to create a groundtruth file.\n'); return; end
disp('Using manually curated clusters');

if isempty(Sclu), Sclu = Sgt; disp('Self-check, expecting 100% accuracy'); end
fprintf('verifying cluster...\n'); 

[mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = clusterVerify(Sgt.viClu, Sgt.viTime, Sclu.viClu, Sclu.viTime);  %Sgt.viTime
% nClu_gt = Sgt.
% for iClu=:nClu_gt


end %func


function batch_(vcFile_batch, vcCommand)
% batch process parameter files

if nargin<2, vcCommand=[]; end
if isempty(vcCommand), vcCommand='spikesort'; end
csFiles_prm = importdata(vcFile_batch);
for iFile=1:numel(csFiles_prm)
    try
        vcFile_prm1 = csFiles_prm{iFile};
        jrc('clear');
        jrclust(vcCommand, vcFile_prm1);
        if isempty(strfind(vcCommand, 'verify'))
            jrclust('verify', vcFile_prm1); 
        end
    catch
        disp(lasterr());
    end
end %for

end %func


function import_abf_(vcFile_abf, thresh1)
% Generate GT file and .bin and .meta
fprintf('Converting from ABF to WHISPER format (.bin, .meta and _gt.mat)\n\t'); t1=tic;
[mrWav, sampling_interval_us, Sfile]=abfload(vcFile_abf);
% P.nChans = size(mnWav,2);
P.nChans = 3;
P.uV_per_bit = 6.103515761424;
P.sRateHz = 1e6 / double(sampling_interval_us);
mrWav(1:20000,:) = []; %cut off first 1 sec

fprintf('\n\ttook %0.1fs.\n', toc(t1));
vrInt = mrWav(:,3);
vrInt = diff_filter(vrInt);
mrWav = int16(mrWav(:, [1,2,4]) * 1000 / P.uV_per_bit); 
% output bin
vcFile_bin = subsFileExt(vcFile_abf, '.bin');
fid = fopen(vcFile_bin, 'W');
fwrite(fid, mrWav', 'int16');
fclose(fid);

int_med = median(vrInt);
if isempty(thresh1)
    thresh1 = int_med + (max(vrInt)-int_med) * .5; 
else
    thresh1 = str2double(thresh1);
end
[Sgt.viTime, ~] = findPeak(vrInt, thresh1);
findPeak(vrInt, thresh1); figure; plot(mrWav);
Sgt.viClu = ones(size(Sgt.viTime));
vcFile_gt = subsFileExt(vcFile_abf, '_gt.mat');
save(vcFile_gt, '-struct', 'Sgt');

% write .meta file
% vcFile_meta = subsFileExt(vcFile, '.meta');
% fid = fopen(vcFile_meta, 'W');
% fprintf(fid, 'niMNGain=1\n');
% fprintf(fid, 'niSampRate=%f\n', P.sRateHz); %intan hardware default. in Smeta.header
% fprintf(fid, 'niAiRangeMax=.5\n'); %intan hardware default. in Smeta.header
% fprintf(fid, 'niAiRangeMin=-.5\n'); %intan hardware default. in Smeta.header
% fprintf(fid, 'nSavedChans=%d\n', P.nChans); %intan hardware default. in Smeta.header
% fclose(fid);
% write meta file and bin file

fprintf('Exported to %s, %s\n', vcFile_bin, vcFile_gt);
end


function import_ncs_(vcFile_ncs)
% import neuralynx file
[mrWav, Smeta, vcFile_bin] = ncs2dat(vcFile_ncs);
% save Smeta to .meta file
fprintf('Exported to %s\n.', vcFile_bin);
disp(Smeta)
end %func


function S_metric = quality_metric_(arg1, arg2, arg3)
% quality accessment of the sorted unit
% S_metric = quality_metric_(vcFile_prm)
% S_metric = quality_metric_(Sevt, Sclu, mrWav)
% no output: write to a file

global mrWav
if nargin==1
    vcFile_prm = arg1;
    P = loadParam_(vcFile_prm);
    % event metric
    Sclu = load(subsFileExt(vcFile_prm, '_clu.mat'));
    Sevt = load(subsFileExt(vcFile_prm, '_evt.mat'));
    mrWav = load_spk_(vcFile_prm);
elseif nargin==3
    [Sevt, Sclu, mrWav] = deal(arg1, arg2, arg3);
    P = Sclu.P; 
    vcFile_prm = P.vcFile_prm;
end
t_dur = single(Sclu.viTime(end)) / (Sclu.P.sRateHz); % max(cellfun(@(x)x(end), Sclu11.cviSpk_site, 'UniformOutput', 1)) / Sevt11.P.sRateHz;            
vrVp_evt = abs(single(Sevt.vrSpk) * P.uV_per_bit); %negative peak
vrVrms_site = Sevt.vrThresh_uV / P.qqFactor;
vrSnr_evt = vrVp_evt ./ vrVrms_site(Sevt.viSite);
vrRate_site = hist(Sevt.viSite, 1:numel(Sevt.P.viSite2Chan))' / t_dur;

% cluster metric
P.nDiff_filt=0; %do not integrate the waveform
P.spkLim = P.spkLim_fet;
trWav_clu = Sclu_meanWav(Sclu, mrWav, P); %compare before and after integration?
mrVmax_clu = shiftdim(max(trWav_clu,[],1));
mrVmin_clu = shiftdim(min(trWav_clu,[],1));
mrVpp_clu = mrVmax_clu - mrVmin_clu;
[vrVp_clu, viSite_clu] = min(mrVmin_clu, [], 1); 
vrVp_clu = abs(vrVp_clu(:));
vrVpp_clu = get_mr(mrVpp_clu, viSite_clu, []);
vrSnr_clu = vrVp_clu ./ vrVrms_site(viSite_clu); 
vrSnr_pp_clu = vrVpp_clu ./ vrVrms_site(viSite_clu); 
vnSpk_clu = Sclu.vnSpk_clu;
vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -Sevt.vrThresh_uV),1)';

% write to file
S_metric = makeStruct(vrSnr_evt, vrVrms_site, vrRate_site, ...
    viSite_clu, vrSnr_clu, vrSnr_pp_clu, t_dur, trWav_clu, vnSpk_clu, vnSite_clu, P);

if nargout==0
    vcFile_metric = subsFileExt(vcFile_prm, '_metric.mat');
    save(vcFile_metric, '-struct', 'S_metric');
    fprintf('\n[%s] Saved to %s\n', datestr(now), vcFile_metric);
    assignWorkspace(S_metric);
end
end %func


function feature_map_(vcFile_prm)
% plot features as a function of spike position
Sevt = load_evt(vcFile_prm);
Sclu = load_clu_(vcFile_prm);
vrY_spk = squeeze(Sevt.trFet(1,1,:));
mrA_spk = squeeze(Sevt.trFet(1,3:end,:));

% func1 = @(x)sqrt(abs(x));
func1 = @(x)x;
vrA_left_spk = func1(mrA_spk(3,:));
vrA_right_spk = func1(mrA_spk(8,:));

viPlot = find(vrY_spk >= 800 & vrY_spk < 1000);
viPlot = viPlot(1:4:end);
func2 = @(x)x-min(x)+1;
viClu1 = func2(Sclu.viClu(viPlot));
clu_map = randperm(max(viClu1));
viClu1 = clu_map(viClu1);
if 0
    figure; 
    AX(1)=subplot(221); scatter(vrY_spk(viPlot), vrA_left_spk(viPlot), 3, viClu1, 'filled'); grid on; xlabel('Spike position (um)'); ylabel('Sqrt Ampl');
    hold on; hPlot=plot(vrY_spk(Sclu.icl), vrA_left_spk(Sclu.icl), 'ro', 'MarkerSize', 20);
    AX(2)=subplot(222); scatter(vrY_spk(viPlot), vrA_right_spk(viPlot), 3, viClu1, 'filled'); grid on; xlabel('Spike position (um)'); ylabel('Sqrt Ampl');
    hold on; hPlot=plot(vrY_spk(Sclu.icl), vrA_right_spk(Sclu.icl), 'ro', 'MarkerSize', 20);
    linkaxes(AX,'xy');

    xlim([800 1000]);
else
    figure; 
    AX(1)=subplot(2,2,1); plot(vrY_spk(viPlot), vrA_left_spk(viPlot), 'k.', 'MarkerSize', 5); grid on; xlabel('Spike position (um)'); ylabel('amp ch3');
%     hold on; hPlot=plot(vrY_spk(Sclu.icl), vrA_left_spk(Sclu.icl), 'ro', 'MarkerSize', 20);
    AX(2)=subplot(2,2,3); plot(vrY_spk(viPlot), vrA_right_spk(viPlot), 'k.', 'MarkerSize', 5); grid on; xlabel('Spike position (um)'); ylabel('amp ch7');
%     hold on; hPlot=plot(vrY_spk(Sclu.icl), vrA_right_spk(Sclu.icl), 'ro', 'MarkerSize', 20);   
    linkaxes(AX,'xy');
    axis([800 1000 0 1500]);
    subplot(2,2,[2,4]); plot(vrA_left_spk(viPlot), vrA_right_spk(viPlot), 'k.', 'MarkerSize', 5); grid on; xlabel('amp ch3'); ylabel('amp ch7');
    axis([0 1500 0 1500]);
end
colormap jet
end %func


function import_celldb_(vcFile_csv)

%read file
fid=fopen(vcFile_csv); 
C =textscan(fid, '%d %s %f %f %f %f %f %s %s', 'Delimiter', ',' ,'HeaderLines', 1); 
fclose(fid);

csType_cell = cellfun(@(x)x(2:end-1), C{2}, 'UniformOutput', 0);
mrPosAng_cell = cell2mat(C([3:end-2]));
viId_cell = C{1};

vcFile_out = strrep(vcFile_csv, '.csv', '.mat');
save(vcFile_out, 'csType_cell', 'mrPosAng_cell', 'viId_cell');
fprintf('Wrote to %s\n', vcFile_out);
end


function export_imec_sync_(vcFiles_prm)
% sync output for IMEC

try
    P = file2struct(vcFiles_prm);
    vcFile_bin = P.vcFile;    
    fid = fopen(vcFile_bin, 'r');
    vnSync = fread(fid, inf, '*uint16'); fclose(fid);
    vnSync = vnSync(P.nChans:P.nChans:end); %subsample sync channel

    % @TODO: skip more efficient memorywise
    % fseek(fid, (P.nChans-1)*2, 'bof');
    % vnSync = fread(fid, inf, '*uint16', 2*(P.nChans-1)); fclose(fid);

    assignWorkspace(vnSync);
catch
    fprintf(2, 'error exporting sync: %s\n', lasterr());
end
% viSync = find(diff(vnSync)~=0)+1;
% .nev file


end


%--------------------------------------------------------------------------
function batch_verify_(vcFile_batch, vcCommand)
% batch process parameter files (.batch) file
% Example
%   jrclust batch-verify skip my.batch
%       just does the verification plot for all files in .batch file

if nargin<2, vcCommand=[]; end
if isempty(vcCommand), vcCommand='spikesort'; end
csFiles_prm = importdata(vcFile_batch);
if ~strcmpi(vcCommand, 'skip')
    for iFile=1:numel(csFiles_prm)
        try
            vcFile_prm1 = csFiles_prm{iFile};
            jrclust('clear');
            jrclust(vcCommand, vcFile_prm1);
            if isempty(strfind(vcCommand, 'verify'))
                jrclust('verify', vcFile_prm1); % try silent verify and collect result
            end
        catch
            disp(lasterr());
        end
    end %for
end

% Collect data
[csSnr, csFp, csFn, cvnSite] = deal(cell(size(csFiles_prm)));
for iFile=1:numel(csFiles_prm)    
    S_score1 = load(strrep(csFiles_prm{iFile}, '.prm', '_score.mat'));    
    csSnr{iFile} = gather(S_score1.vrSnr_min_gt');    
    csFp{iFile} = S_score1.S_score_clu.vrFp;
    csFn{iFile} = S_score1.S_score_clu.vrMiss;
    cvnSite{iFile} = S_score1.vnSite_gt;
%     csNspk{iFile} = arrayfun(@(x)sum(S_score1.Sgt.viClu==x), 1:max(S_score1.Sgt.viClu));   
%     Sclu1 = load(subsFileExt(csFiles_prm{iFile}, '_clu.mat'));
%     csSnr_clu{iFile} = Sclu1.vrSnr_Vmin_clu;
end

[vrSnr, vrFp, vrFn, vnSite] = multifun(@(x)cell2mat(x'), csSnr, csFp, csFn, cvnSite);
disp_score_(vrSnr, vrFp, vrFn, 8);

% Plot
vpp_bin = 2; %25 for vpp

figure;  ax=[];
ax(1)=subplot(311); 
boxplot_jjj(vrFp(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Positive'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
ylim([0, .2]); xlim([0 40]);

ax(2)=subplot(312); 
boxplot_jjj(vrFn(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Negative'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
ylim([0, .2]); xlim([0 40]);
set(gcf,'Color','w');

ax(3)=subplot(313); 
boxplot_jjj(vnSite(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('#sites>thresh'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
linkaxes(ax,'x');
set(gcf,'Color','w');
ylim([0, 16]); 
end %func


%--------------------------------------------------------------------------
function disp_score_(vrSnr, vrFp, vrFn, snr_thresh)
% disp_score_(vrSnr, vrFp, vrFn, snr_thresh)
fprintf('SNR>%d Groundtruth Units\n', snr_thresh);
viSnr = find(vrSnr > snr_thresh);
fprintf('\tfalse positive (%%): '); disp_stats(vrFp(viSnr)*100);
fprintf('\tfalse negative (%%): '); disp_stats(vrFn(viSnr)*100);
end %func


%--------------------------------------------------------------------------
function nBytes = mem_max_(P)
if P.fGpu
    S = gpuDevice(1);
    nBytes = floor(S.AvailableMemory());
else
    S = memory();
    nBytes = floor(S.MaxPossibleArrayBytes());
end
end %func