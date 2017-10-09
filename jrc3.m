%--------------------------------------------------------------------------
% JRCLUST v3
% James Jun
 
function varargout = jrc3(vcCmd, vcArg1, vcArg2, vcArg3, vcArg4, vcArg5)
% Memory-efficient version. Rewrote from scratch and minimalistic
% P is static and loaded from file. For dynamic property set S0
persistent vcFile_prm_
% global tnWav_spk viSite_spk viTime_spk_ P

% input parse
if nargin<2, vcArg1=''; end
if nargin<3, vcArg2=''; end
if nargin<4, vcArg3=''; end
if nargin<5, vcArg4=''; end
if nargin<6, vcArg5=''; end
if nargin==0, vcCmd = 'help'; end
setpath_();

% Basic commands
fExit = 1;
switch lower(vcCmd)
    % No arguments
    case 'version', jrc_version_(vcArg1);
    case {'help', '-h', '?', '--help'}, help_(vcArg1); about_();
    case 'about', about_();
    case 'clear', clear_(vcArg1);
    case 'doc', doc_('JRCLUST manual.pdf');
    case 'doc-edit', doc_('JRCLUST manual.docx');
    case 'update', update_(vcArg1);
    case 'git-pull', git_pull_(vcArg1);
    case 'install', install_();
    case 'commit', commit_(vcArg1);
    case 'wiki', wiki_(vcArg1);
    case 'wiki-download', wiki_download_();
    case 'gui', gui_(vcArg1, vcFile_prm_);
    case 'issue', issue_('post');
        
    case 'which', return;    
    case 'download', download_(vcArg1);
    case {'makeprm', 'createprm'}
        vcFile_prm_ = makeprm_(vcArg1, vcArg2, 1, vcArg3);
        if nargout>0, varargout{1} = vcFile_prm_; end
    case 'makeprm-f', makeprm_(vcArg1, vcArg2, 0, vcArg3);
    case 'import-tsf', import_tsf_(vcArg1);
    case 'import-jrc1', import_jrc1_(vcArg1);
    case 'export-jrc1', export_jrc1_(vcArg1);
    case 'import-intan', vcFile_prm_ = import_intan_(vcArg1, vcArg2, vcArg3); return;
    case {'import-nsx', 'import-ns5'}, vcFile_prm_ = import_nsx_(vcArg1, vcArg2, vcArg3); return;
    case 'nsx-info', [~, ~, S_file] = nsx_info_(vcArg1); assignWorkspace_(S_file); return;
    case 'load-nsx', load_nsx_(vcArg1); return;
    case 'load-bin'
        mnWav = load_bin_(vcArg1, vcArg2); 
        assignWorkspace_(mnWav);
    case 'import-gt', import_gt_silico_(vcArg1);   
    case 'unit-test', unit_test_(vcArg1, vcArg2, vcArg3);    
    case 'compile', compile_cuda_(); 
    case 'compile-ksort', compile_ksort_();
    case 'test', varargout{1} = test_(vcArg1, vcArg2, vcArg3, vcArg4);
    case 'export', export_(vcArg1, vcArg2, vcArg3);    
    case {'dependencies', 'toolbox', 'toolboxes'}, disp_dependencies_();
    otherwise, fExit = 0;
end
if fExit, return; end


% Commands requiring .prm file
if nargin>=2
    vcFile_prm = vcArg1; 
    vcFile_prm_ = vcFile_prm; 
else
    vcFile_prm = vcFile_prm_;    
end
if isempty(vcFile_prm), disp('Please specify .prm file.'); return; end
if isempty(vcArg1) && ~isempty(vcFile_prm), disp(['Working on ', vcFile_prm]); end
fExit = 1;
switch lower(vcCmd)    
    case 'probe', probe_(vcFile_prm);
    case 'edit', edit_(vcFile_prm); 
    case 'batch', batch_(vcArg1, vcArg2); 
    case 'batch-mat', batch_mat_(vcArg1, vcArg2); %text file containing binary files and template file
    case {'batch-verify', 'batch-validate'}, batch_verify_(vcArg1, vcArg2); 
    case {'batch-plot', 'batch-activity'}, batch_plot_(vcArg1, vcArg2); 
    case 'describe', describe_(vcFile_prm); 
    case 'import-silico', import_silico_(vcFile_prm, 0); 
    case 'import-silico-sort', import_silico_(vcFile_prm, 1); 
    case {'import-kilosort', 'import-ksort'}, import_ksort_(vcFile_prm, 0); 
    case {'import-kilosort-sort', 'import-ksort-sort'}, import_ksort_(vcFile_prm, 1);  
    case {'kilosort', 'ksort'}, kilosort_(vcFile_prm); import_ksort_(vcFile_prm, 0); 
    case 'export-imec-sync', export_imec_sync_(vcFile_prm);
    case 'export-prm', export_prm_(vcFile_prm, vcArg2);        
    case 'dir', 
        if any(vcFile_prm=='*')
            dir_files_(vcFile_prm, vcArg2, vcArg3);
        else
            fExit = 0;
        end
    otherwise
        fExit = 0;
end
if fExit, return; end

% Reset GPU
% try gpuDevice(1); disp('GPU device reset'); catch, end; %optional

% Commands requiring P structure
if ~matchFileExt_(vcFile_prm, '.prm'), fprintf(2, 'Must provide .prm file\n'); return ;end
P = loadParam_(vcFile_prm); 
if isempty(P), return; end    
fError = 0;
switch lower(vcCmd)    
    case 'preview', preview_(P); 
    case 'preview-test', preview_(P, 1); gui_test_(P, 'Fig_preview');
    case 'traces', traces_(P, 0, vcArg2);
    case 'dir', dir_files_(P.csFile_merge);
    case 'traces-test'
        traces_(P, 1); traces_test_(P);           
    case {'full', 'all'}
        detect_(P); sort_(P); describe_(P.vcFile_prm); manual_(P); return;      
    case {'spikesort', 'detectsort', 'detect-sort', 'spikesort-verify', 'spikesort-validate', 'spikesort-manual', 'detectsort-manual'}
        detect_(P); sort_(P); describe_(P.vcFile_prm);
    case {'detect', 'spikedetect'}
        detect_(P); describe_(P.vcFile_prm);
    case {'sort', 'cluster', 'clust', 'sort-verify', 'sort-validate', 'sort-manual'}
        if ~is_detected_(P), detect_(P); end
        sort_(P); describe_(P.vcFile_prm);
    case {'auto', 'auto-verify', 'auto-manual'}
        auto_(P); describe_(P.vcFile_prm);
    case 'manual-test'
        manual_(P, 'debug'); manual_test_(P); return;          
    case 'manual-test-menu'
        manual_(P, 'debug'); manual_test_(P, 'Menu'); return;                   
    case {'kilosort-verify', 'ksort-verify'}
        kilosort_(P); import_ksort_(P); describe_(P.vcFile_prm); 
    case {'export-wav', 'wav'} % load raw and assign workspace
        mnWav = load_file_(P.vcFile, [], P);
        assignWorkspace_(mnWav);
    case 'export-spk'
        S0 = get0_();
        trSpkWav = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
        assignWorkspace_(trSpkWav);        
    case 'export-raw'
        S0 = get0_();
        trWav_raw = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_spk);
        assignWorkspace_(trWav_raw);  
    case {'export-spkwav', 'spkwav'}, export_spkwav_(P, vcArg2); % export spike waveforms
    case {'export-spkwav-diff', 'spkwav-diff'}, export_spkwav_(P, vcArg2, 1); % export spike waveforms        
    case 'export-spkamp', export_spkamp_(P, vcArg2); %export microvolt unit
    case {'export-csv', 'exportcsv'}, export_csv_(P);
    case {'activity', 'plot-activity'}, plot_activity_(P);    
    case {'export-fet', 'export-features', 'export-feature'}, export_fet_(P);
    case 'export-diff', export_diff_(P); %spatial differentiation for two column probe
    case 'drift', plot_drift_(P);  
    case 'plot-rd', plot_rd_(P);
    otherwise, fError = 1;
end %switch

% supports compound commands (ie. 'sort-verify', 'sort-manual').
if contains_(lower(vcCmd), {'verify', 'validate'})
    if ~is_detected_(P), detect_(P); end
    if ~is_sorted_(P), sort_(P); end    
    validate_(P);
elseif contains_(lower(vcCmd), {'manual',' gui', 'ui'})
    manual_(P);
elseif fError
    help_();
end
end %func


%--------------------------------------------------------------------------
function flag = contains_(vc, cs)
% check if vs contains any of cs
if ischar(cs), cs = {cs}; end
flag = 0;
for i = 1:numel(cs)
    if ~isempty(strfind(vc, cs{i}))
        flag = 1; 
        break; 
    end
end

end %func


%--------------------------------------------------------------------------
function [viSite_spk2, viSite_spk3] = find_site_spk23_(tnWav_spk, viSite_spk, P)
% find second min, excl local ref sites
fUse_min = 0;
imin0 = 1 - P.spkLim(1);
viSites2 = 2:(2*P.maxSite+1-P.nSites_ref);
miSites2 = P.miSites(viSites2,viSite_spk);
%[~,viSite_spk2] = min(squeeze_(tnWav_spk(imin0,viSites2,:)));
tnWav_spk2 = tnWav_spk(:,viSites2,:);
if fUse_min
    mnMin_spk = squeeze_(min(tnWav_spk2));
else
%     mnMin_spk = -squeeze_(std(single(tnWav_spk(:,viSites2,:))));
    mnMin_spk = squeeze_(min(tnWav_spk2) - max(tnWav_spk2)); % use Vpp to determine second peak site
end
if nargout==1
    [~, viSite_spk] = min(mnMin_spk);
    viSite_spk2 = int32(mr2vr_sub2ind_(miSites2, viSite_spk, []));
else
    [~, miSite_spk2] = sort(mnMin_spk, 'ascend');
    viSite_spk2 = int32(mr2vr_sub2ind_(miSites2, miSite_spk2(1,:), []));
    viSite_spk3 = int32(mr2vr_sub2ind_(miSites2, miSite_spk2(2,:), []));
end
end %func


%--------------------------------------------------------------------------
function tnWav_spk1 = mn2tn_wav_spk2_(mnWav1, viSite_spk, viTime_spk, P)
nSpks = numel(viSite_spk);
nSites = numel(P.viSite2Chan);
spkLim_wav = P.spkLim;
nSites_spk = (P.maxSite * 2) + 1;
tnWav_spk1 = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'int16');
for iSite = 1:nSites
    viiSpk11 = find(viSite_spk == iSite);
    if isempty(viiSpk11), continue; end
    viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
    viSite11 = P.miSites(:,iSite);
    tnWav_spk1(:,:,viiSpk11) = permute(gather_(mr2tr3_(mnWav1, spkLim_wav, viTime_spk11, viSite11)), [1,3,2]); %raw
end
end %func


%--------------------------------------------------------------------------
function csHelp = help_(vcCommand)
if nargin<1, vcCommand = ''; end
if ~isempty(vcCommand), wiki_(vcCommand); return; end

csHelp = {...
    ''; 
    'Usage: jrc command arg1 arg2 ...';
    '  jrc calls the latest jrclust version (currently jrc3.m)';
    '  You can also run "jrc3 command ..." to explicitly specify the version number.';
    '';
    '[Documentation and help]';
    '  jrc help';
    '    Display a help menu'; 
    '  jrc doc';
    '    Open a help document (pdf)';         
    '  jrc version';
    '    Display the version number and the updated date';            
    '  jrc about';
    '    Display about information';       
    '  jrc wiki';
    '    Open a JRCLUST Wiki on GitHub';     
    '  jrc wiki-download';
    '    Download the JRCLUST Wiki on GitHub to ./wiki/';    
    '  jrc issue';
    '    Post an issue at GitHub (log-in with your GitHub account)';        
    '';
    '[Main commands]';
    '  jrc edit (myparam.prm)';
    '    Edit .prm file currently working on'; 
    '  jrc clear';
    '    Clear cache';
    '  jrc clear myparam.prm';
    '    Delete previous results (files: _jrc.mat, _spkwav.jrc, _spkraw.jrc, _spkfet.jrc)';        
    '  jrc spikesort myparams.prm';
    '    Run the whole suite (spike detection and clustering) ';
    '  jrc detect myparams.prm';
    '    Run spike detection and extract spike waveforms';
    '    Output files: _jrc.mat, _spkwav.jrc (filtered spike waveforms), _spkraw.jrc (raw spike waveforms), _spkfet.jrc (features)';
    '  jrc sort myparams.prm';
    '    Cluster spikes (after spike detection)';
    '    Output files: _jrc.mat';
    '  jrc auto myparams.prm';
    '    Recluster spikes after updating post-clustering paramters';
    '  jrc download sample';
    '    Download sample data from Neuropix phase 2 probe';
    '  jrc probe {myprobe.prb, myparams.prm}';
    '    Plot probe layout'; 
    '  jrc makeprm myrecording.bin myprobe.prb';
    '    create a new parameter file based on the default template file (default.prm) and probe file';
    '  jrc makeprm myrecording.bin myprobe.prb mytemplate.prm';
    '    create a new parameter file based on the specified template file and probe file';    
    '  jrc traces myparams.prm';
    '    Displays raw trace';
    '  jrc describe myparams.prm';
    '    Display information about a clu dataset ';
    '  jrc manual myparams.prm';
    '    Run the manual clustering GUI ';
    '  jrc auto-manual myparams.prm';
    '    Run the auto clustering and do the manual clustering next';        
    '  jrc plot-activity myparams.prm';
    '    Show firing rate as a function of time and depth';         
    '  jrc verify myparams.prm';
    '    Compares against ground truth file (_gt.mat)';
    '  jrc drift myparams.prm';
    '    Visualize drift.';       
    '  jrc plot-rd myparams.prm';
    '    Rho vs. Delta plot used in DPCLUS clustering (Rodriguez-Laio).';      
    '  jrc batch myparam.batch (command)';
    '    Batch process list of prm files';         
    '  jrc batch myparam.batch template.prm';
    '    Batch process list of .bin files using a template prm file';             
    '';
    '[Import and export]';
    '  jrc export myparams.prm';
    '    Export the global struct (S0) to the workspace. This is also contained in _jrc.mat output file.';    
    '  jrc export-csv myparams.prm';
    '    Export clustered information to a csv file (spike time, cluster #, max site#)';
    '  jrc export-jrc1 myparams.prm';
    '    Export to version 1 format (write to _evt.mat and _clu.mat)';
    '  jrc import-jrc1 myparams.prm';
    '    Import from version 1 format';    
    '  jrc export-imec-sync myparams.prm';
    '    Export Sync channel (uint16) to the workspace (vnSync)';            
    '  jrc export-wav myparams.prm';
    '    Export the entire raw traces (mnWav) to the Workpace';
    '  jrc export-spkwav myparams.prm (clu#)';
    '    Export spike waveforms organized by clusters to the Workspace';
    '  jrc export-spk myparams.prm';
    '    Export spike waveforms to the Workspace, sorted by the time of spike';    
    '  jrc export-spkamp myparams.prm (clu#)';
    '    Export spike amplitudes to the Workspace'; 
    '  jrc export-fet myparams.prm';
    '    Export feature matrix (trFet) and sites (miFet_sites) to the Workspace';     
    '  jrc export-prm myparams.prm (myparam_full.prm)';
    '    Export complete list of parameters to a new file.';
    '    If the second argument is omitted, the current parameter file is updated.';
    '    If some parameter values are missing, they will be copied from the default template file (default.prm).';
    '  jrc import-intan intanrec*.dat myprobe.prb';
    '    Import from intan recordings, which saved each channel as -A###.dat file in int16 format.';        
    '    Combine .dat files to a single .bin file and saves to the directory above.';
    '    Generate .prm file by combining the .bin and .prb file names (e.g. binfile_prbfile.prm).';    
    '  jrc import-nsx myrec.ns5 myprobe.prb';
    '    Imports Neuroshare format and export the analog channels to .bin file.';
    '    Generates .prm file by combining the .bin and .prb file names (e.g. binfile_prbfile.prm).';
    '';
    '[Sorting multiple recordings together]';
    '  jrc dir myparam.prm'; 
    '    List all recording files to be clustered together (csFile_merge)';
    '  jrc traces myparam.prm';
    '    List all recording files and select which one to display';
    '  jrc traces myparam.prm File#';
    '    Direcly specify the file number to display';
    '';
    '[Deployment]';
    '  jrc update';
    '    Update code by copying from the dropbox location (specified in user.cfg or default.cfg)';    
    '  jrc unit-test';
    '    Run a suite of unit teste.';       
    '  jrc install';
    '    Install jrc by compiling codes';    
    '  jrc compile';
    '    Recompile CUDA code (GPU codes, *.cu)';     
    '';
};
if nargout==0, disp_cs_(csHelp); end
end %func


%--------------------------------------------------------------------------
function disp_cs_(cs)
% display cell string
cellfun(@(s)fprintf('%s\n',s), cs);
end %func


%--------------------------------------------------------------------------
function update_(vcFile)
% update_(vcFile) %update specific files
% update_() %update files in default.cfg\sync_list
fCompile_ksort = 0;

S_cfg = read_cfg_();
vcSource = S_cfg.path_dropbox;
sprintf('copyfile ''%s\\%s'' .\\ f;', vcSource, 'default.cfg'); % copy default
pause(.5); S_cfg = read_cfg_(); % reload default
if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
vcBackup = S_cfg.path_backup;
mkdir_(vcBackup);

t1 = tic;
fprintf('Updating from %s...\n', vcSource);
csCopy = ifeq_(isempty(vcFile), S_cfg.sync_list, {vcFile});   
for iCopy = 1:numel(csCopy)
    try_eval_(sprintf('copyfile ''%s\\%s'' .\\ f;', vcSource, csCopy{iCopy}), 0);
    try_eval_(sprintf('copyfile ''.\\%s'' ''%s\\'' f;', csCopy{iCopy}, vcBackup), 0);
end

% Compile CUDA code
fCompile = isempty(vcFile);
try
    if fCompile, compile_cuda_(S_cfg); end
catch
    fprintf(2, 'CUDA code compilation error.\n');
end

% Copy kilosort and compile code
if fCompile_ksort
    try 
        mkdir_('./kilosort');
        try_eval_(sprintf('copyfile ''%s\\kilosort\\*'' .\\kilosort\\ f;', vcSource), 0);
        if fCompile, compile_ksort_(); end
    catch
        disperr_();
    end
end

fprintf('Updated, took %0.1fs.', toc(t1));
fprintf('\tPrevious files backed up to %s\n', vcBackup);
edit change_log.txt
end %func


%--------------------------------------------------------------------------
function probe_(vcFile_prb)
% if nargin<1, vcFile_prb='imec2.prb'; end

% set prb file
if matchFileExt_(vcFile_prb, {'.bin', '.dat'})
    vcFile_prb = subsFileExt_(vcFile_prb, '.prm');
end
if matchFileExt_(vcFile_prb, '.prm')
    vcFile_prm = vcFile_prb;
    P = loadParam_(vcFile_prm);
    vcFile_prb = P.probe_file;
    if ~exist(vcFile_prb, 'file')
        vcFile_prb = replacePath_(vcFile_prb, vcFile_prm);
    end
end
vcFile_prb = find_prb_(vcFile_prb); 
S_prb = file2struct_(vcFile_prb);
if ~isfield(S_prb, 'shank'), S_prb.shank = ones(size(S_prb.channels)); end

% hFig = figure; hold on;
hFig = create_figure_('FigProbe', [0 0 .5 1], vcFile_prb);
% viShank = unique(S_prb.shank);
% vcColor_shank = 'kbgrcm'; % up to 
% for iShank=1:numel(viShank)
%     viSite1 = find(S_prb.shank == iShank);
hPatch = plot_probe_(S_prb.geometry, S_prb.pad, S_prb.channels, S_prb.shank);
%     iShank1 = mod(iShank-1, numel(vcColor_shank))+1;
%     set(hPatch, 'EdgeColor', vcColor_shank(iShank1));
% end
% plot_probe_(mrSiteXY, vrSiteHW, viSite2Chan, vrVpp, hFig)
% vrPos0 = get(0, 'ScreenSize'); 
% set(hFig, 'OuterPosition', vrPos0, 'Color', 'w');
axis equal;
% set(hFig, 'Name', vcFile_prb, 'NumberTitle', 'off');
edit(vcFile_prb); %show probe file
figure(hFig);
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Find
function vcFile_prb = find_prb_(vcFile_prb)
% Find a prb file
vcFile_prb = search_file_(vcFile_prb, [jrcpath_(), 'prb', filesep()]);
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_full = jrcpath_(vcFile, fConditional)
% make it a jrc path
% Add a full path if the file doesn't exist in the current folder
% 

if nargin<1, vcFile = ''; end
if nargin<2, fConditional = 0; end

jrcpath = fileparts(mfilename('fullpath'));
if fConditional
    if exist(vcFile, 'file') == 2
        vcFile_full = vcFile;
        return;
    end
end
vcFile_full = [jrcpath, filesep(), vcFile];
% if exist(vcFile_full, 'file') ~= 2
%     vcFile_full = [];
% end
end %func


%--------------------------------------------------------------------------
function download_(vcMode)
S_cfg = read_cfg_();
% if strcmpi(pwd(), S_cfg.path_alpha), disp('cannot overwrite alpha'); return; end
switch lower(vcMode)   
    case {'sample', 'neuropix2', 'neuropixels2', 'phase2', 'phaseii'}
        csLink = S_cfg.path_sample_phase2;
    case {'sample3', 'neuropix3' 'neuropixels3', 'phase3', 'phaseiii'}
        csLink = S_cfg.path_sample_phase3;
    otherwise
        disp('Invalid selection. Try "jrc download sample or jrc download sample3".');
        return;
end %switch

t1 = tic;
fprintf('Downloading sample files. This can take up to several minutes.\n');
vlSuccess = download_files_(csLink);
fprintf('\t%d/%d files downloaded. Took %0.1fs\n', ...
    sum(vlSuccess), numel(vlSuccess), toc(t1));
end %func


%--------------------------------------------------------------------------
function keyPressFcn_Fig_traces_(hFig, event)
% 2017/6/22 James Jun: Added nTime_traces multiview

global mnWav1 mrWav1 mnWav
[S0, P] = get0_();
S_fig = get(hFig, 'UserData');
factor = 1 + 3 * key_modifier_(event, 'shift');
nSites = numel(P.viSite2Chan);

switch lower(event.Key)
    case 'h', msgbox_(S_fig.csHelp, 1);
        
    case {'uparrow', 'downarrow'}
        if isfield(S_fig, 'chSpk')
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot, S_fig.chSpk);
        else
            S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, S_fig.hPlot);
        end
        title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));        
        set(hFig, 'UserData', S_fig);
        
    case {'leftarrow', 'rightarrow', 'j', 'home', 'end'}
        switch lower(event.Key)
            case 'leftarrow'
                nlim_bin = S_fig.nlim_bin - (S_fig.nLoad_bin) * factor; %no overlap
                if nlim_bin(1)<1
                    msgbox_('Beginning of file', 1); 
                    nlim_bin = [1, S_fig.nLoad_bin]; 
                end
            case 'rightarrow'
                nlim_bin = S_fig.nlim_bin + (S_fig.nLoad_bin + 1) * factor; %no overlap
                if nlim_bin(2) > S_fig.nSamples_bin
                    msgbox_('End of file', 1); 
                    nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
                end
            case 'home' %beginning of file
                nlim_bin = [1, S_fig.nLoad_bin];
            case 'end' %end of file
                nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
            case 'j'
                vcAns = inputdlg_('Go to time (s)', 'Jump to time', 1, {'0'});
                if isempty(vcAns), return; end
                try
                    nlim_bin = round(str2double(vcAns)*P.sRateHz) + [1, S_fig.nLoad_bin];
                catch
                    return;
                end
        end %switch
        nTime_traces = get_(P, 'nTime_traces');
        [cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, S_fig.nSamples_bin, nTime_traces); 
        if P.fTranspose_bin
            fseek_(S_fig.fid_bin, nlim_bin(1), P);
            if nTime_traces > 1
                mnWav1 = load_bin_multi_(S_fig.fid_bin, cvn_lim_bin, P)';
            else
                mnWav1 = load_bin_(S_fig.fid_bin, P.vcDataType, [P.nChans, S_fig.nLoad_bin])';            
            end
        else
            mnWav1 = mnWav(viRange_bin, :);
        end
        S_fig.nlim_bin = nlim_bin;
        set_fig_(hFig, S_fig);
        Fig_traces_plot_(1); %redraw
        
    case 'f' %apply filter
        S_fig.vcFilter = str_toggle_(S_fig.vcFilter, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();        
        
    case 'g' %grid toggle on/off
        S_fig.vcGrid = str_toggle_(S_fig.vcGrid, 'on', 'off');
        grid(S_fig.hAx, S_fig.vcGrid);
        set(hFig, 'UserData', S_fig);
        
    case 'r' %reset view
        fig_traces_reset_(S_fig);        
        
    case 'e' %export current view        
        assignWorkspace_(mnWav1, mrWav1);
        disp('mnWav1: raw traces, mrWav1: filtered traces');
        
    case 'p' %power spectrum
        iSite_show = inputdlg_num_(sprintf('Site# to show (1-%d, 0 for all)', nSites), 'Site#', 0);
        if isnan(iSite_show), return; end        
        hFig = create_figure_('FigPsd', [.5 0 .5 1], P.vcFile_prm, 1, 1); %show to the right
        % ask user which channels to plot
        if iSite_show>0
            mrWav2 = mrWav1(:, iSite_show);
        else
            mrWav2 = mrWav1;
        end
        plotMedPower_(mrWav2, 'sRateHz', P.sRateHz/P.nSkip_show, 'viChanExcl', P.viSiteZero);
        
    case 's' %show/hide spikes
        S_fig.vcSpikes = str_toggle_(S_fig.vcSpikes, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();
        
    case 't' %show/hide traces
        S_fig.vcTraces = str_toggle_(S_fig.vcTraces, 'on', 'off');
        set_fig_(hFig, S_fig);
        Fig_traces_plot_();         
        
    case 'c' %channel query
        msgbox_('Draw a rectangle', 1);
        hRect = imrect_();
        if isempty(hRect), return ;end
        vrPos_rect = getPosition(hRect);            
        S_plot = get(S_fig.hPlot, 'UserData');
        vrX = get(S_fig.hPlot, 'XData');
        vrY = get(S_fig.hPlot, 'YData');
        viIndex = find(vrX >= vrPos_rect(1) & vrX <= sum(vrPos_rect([1,3])) & vrY >= vrPos_rect(2) & vrY <= sum(vrPos_rect([2,4])));
        if isempty(viIndex), delete_multi_(hRect); return; end
        index_plot = round(median(viIndex));
        [time1, iSite] = ind2sub(size(mrWav1), index_plot);        
        mrX = reshape(vrX, S_plot.dimm);
        mrY = reshape(vrY, S_plot.dimm);
        hold(S_fig.hAx, 'on');
        hPoint = plot(vrX(index_plot), vrY(index_plot), 'r*');
        hLine = plot(S_fig.hAx, mrX(:,iSite), mrY(:,iSite), 'r-');
        hold(S_fig.hAx, 'off');
        iChan = P.viSite2Chan(iSite);
        msgbox_(sprintf('Site: %d/ Chan: %d', iSite, iChan), 1);
        delete_multi_(hRect, hLine, hPoint);
end %return if S_fig didn't change
end %func


%--------------------------------------------------------------------------
function vc = str_toggle_(vc, vc1, vc2)
% toggle vc1 to vc2
if strcmpi(vc, vc1)
    vc = vc2; 
else
    vc = vc1;
end
end %func


%--------------------------------------------------------------------------
function [maxAmp, mrAmp_prev] = change_amp_(event, maxAmp, varargin)
% varargin: plot object to rescale
% Change amplitude scaling 
if nargin<3, hPlot=[]; end
factor = sqrt(2);
if key_modifier_(event, 'shift'), factor = factor ^ 4; end
mrAmp_prev = maxAmp;
if strcmpi(event.Key, 'uparrow')
    maxAmp = maxAmp / factor;
elseif strcmpi(event.Key, 'downarrow')
    maxAmp = maxAmp * factor;
end
for iPlot = 1:numel(varargin)
    try
        multiplot(varargin{iPlot}, maxAmp);
    catch            
    end
end
% handle_fun_(@rescale_plot_, hPlot, maxAmp);
end %func


%--------------------------------------------------------------------------
function flag = key_modifier_(event, vcKey)
% Check for shift, alt, ctrl press
try
    flag = any(strcmpi(event.Modifier, vcKey));
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function val = read_cfg_(vcName, fVerbose)
% read configuration file that stores path to folder
% load from default.cfg but override with user.cfg if it exists
if nargin<2, fVerbose = 1; end

S_cfg = file2struct_('default.cfg');
% end
if exist('user.cfg', 'file')
    S_cfg1 = file2struct_('user.cfg'); %override
    S_cfg = struct_merge_(S_cfg, S_cfg1, {'path_dropbox', 'path_backup', 'default_prm'});
    if fVerbose, fprintf('Configuration loaded from user.cfg.\n'); end
else
    if fVerbose, fprintf('Configuration loaded from default.cfg.\n'); end
end
if nargin==0
    val = S_cfg; 
else
    try
        val = S_cfg.(vcName);
    catch
        disp_err_('read_cfg_');
        switch lower(vcName)
            case 'default_prm'
                val = 'default.prm';
            otherwise
                val = [];
        end
    end
end
end %func


%--------------------------------------------------------------------------
% update jrclust_alpha to jrclust directory
function commit_(vcArg1)
% commit_()
%   Full validation and update
% commit_('log')
%   Update 'change_log.txt' only
% commit_('skip')
%   Skip unit test

if nargin<1, vcArg1=''; end
t1 = tic;
S_cfg = read_cfg_();
% csCopy = S_cfg.sync_list;
% vcDest = S_cfg.path_dropbox;
% vcDest2 = S_cfg.path_dropbox2;
if ~strcmpi(pwd(), S_cfg.path_alpha), disp('must commit from alpha'); return; end

if strcmpi(vcArg1, 'log')
%     sprintf('copyfile change_log.txt ''%s'' f;', S_cfg.path_dropbox);  
    copyfile_('change_log.txt', ...
        struct_get_(S_cfg, ...
            'path_dropbox', 'path_dropbox2', 'path_web', 'path_bitbucket', 'path_github'));
    disp('Commited change_log.txt');
    return;
elseif ~strcmpi(vcArg1, 'skip')
    disp('Running unit tests before commit... ');
    nTests_fail = unit_test_(); %run full unit test
    if nTests_fail > 0
        fprintf(2, 'Commit aborted, %d unit tests failed\n.', nTests_fail);
        return;
    end
else
    fprintf(2, 'Skipping unit test...\n');
end

% commit all
disp(['Commiting to ', S_cfg.path_dropbox]);
delete_files_(find_empty_files_());
copyfile_(S_cfg.sync_list, S_cfg.path_dropbox); %destination root
% copyfile_('./kilosort/*', [S_cfg.path_dropbox, filesep(), 'kilosort', filesep()]); %destination root
copyfile_(S_cfg.csFiles_sample, S_cfg.path_dropbox);

% commit jrc3 related files only
try commit_jrc3_(S_cfg, [S_cfg.path_web, filesep()], 1); catch; disperr_(); end
try commit_jrc3_(S_cfg, [S_cfg.path_bitbucket, filesep()], 0); catch; disperr_(); end
try commit_jrc3_(S_cfg, [S_cfg.path_github, filesep()], 0); catch; disperr_(); end

edit change_log.txt
fprintf('Commited, took %0.1fs.\n', toc(t1));
end


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function cvr = struct_get_(S, varargin)
% Obtain a member of struct
cvr = cell(size(varargin));
for i=1:numel(varargin)
    vcName = varargin{i};
    if isfield(S, vcName)
        cvr{i} = S.(vcName);
    end
end %for
end %func


%--------------------------------------------------------------------------
function detect_(P)
% global tnWav_raw tnWav_spk;
global tnWav_raw tnWav_spk trFet_spk;
runtime_detect = tic;
% Clear memory (S0 is cleared)
set(0, 'UserData', []);

% Detect spikes from files
% [tnWav_raw, tnWav_spk] = deal([]); % clear memory
[tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P);
set(0, 'UserData', S0);

% Save to file
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk); 
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk); 

% measure time
runtime_detect = toc(runtime_detect);
fprintf('Detection took %0.1fs for %s\n', runtime_detect, P.vcFile_prm);
set0_(runtime_detect);

save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
delete_(strrep(P.vcFile_prm, '.prm', '_log.mat')); %delete log file when detecting
end %func


%--------------------------------------------------------------------------
function S0 = sort_(P)
% Extract feature and sort
% global trFet_spk %spike waveform (filtered)
runtime_sort = tic;

S0 = load_cached_(P); % load cached data or from file if exists

% Sort and save
S_clu = fet2clu_(S0, P);
[S_clu, S0] = S_clu_commit_(S_clu, 'sort_');
set0_(P); %, dimm_fet, cvrTime_site, cvrVpp_site, cmrFet_site, P);

% measure time
runtime_sort = toc(runtime_sort);
fprintf('Sorting took %0.1fs for %s\n', runtime_sort, P.vcFile_prm);
set0_(runtime_sort);
S0 = clear_log_(S0);

save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));

% delete the history
end %func


%--------------------------------------------------------------------------
% 8/1/17 JJJ: Bugfix: If S_clu is empty, reload _jrc.mat file
function S0 = load_cached_(P, fLoadWav)
% Load cached data either from RAM or disk
% Usage
% S0 = load_cached_(P)
% S0 = load_cached_(vcFile)
if nargin<2, fLoadWav=1; end

global tnWav_spk tnWav_raw trFet_spk %spike waveform (filtered)
if ischar(P), P = loadParam_(P); end
S0 = get(0, 'UserData'); 
fClear_cache = 1;
if ~isempty(S0)
    if isfield(S0, 'P')
        if strcmpi(P.vcFile_prm, S0.P.vcFile_prm)
            fClear_cache = 0;
        end
    end    
end
if fClear_cache
    S0 = []; tnWav_spk = []; tnWav_raw = []; trFet_spk = []; % clear all
end

% Load from disk
try
    fLoad0 = 0;
    if isempty(S0)
        fLoad0 = 1; 
    elseif isempty(get_(S0, 'S_clu'))
        fLoad0 = 1;
    end
    
    if fLoad0, S0 = load0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat')); end
    S0.P = P;
    if isempty(tnWav_spk) || isempty(tnWav_raw) || isempty(trFet_spk)
        if ~fLoadWav, return; end
        if isempty(S0), return; end %no info
        if isfield(S0, 'dimm_spk')
            tnWav_spk = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), 'int16', S0.dimm_spk);
        end
        if isfield(S0, 'dimm_raw')
            tnWav_raw = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
        end
        if isfield(S0, 'dimm_fet')
            trFet_spk = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
        end        
    end
catch
    disperr_();
end
% S0.P = P; % force S0 to take the latest values from the prm file
% set(0, 'UserData', S0);
end


%--------------------------------------------------------------------------
function varargout = trWav2fet_cov_(trWav2, P)
% tnWav1: nT x nSites_spk x nSpk
% subtract ref
% nSites_spk = 1 + 2 * P.maxSite - P.nSites_ref; % size(tnWav_spk, 2);
vnDelay_fet = get_set_(P, 'vnDelay_fet', [0,3]);

[nT, nSpk, nSites_spk] = size(trWav2);
[cvi2, cvi1] = shift_range_(nT, [], vnDelay_fet);
cmrFet = cell(numel(vnDelay_fet), 1);
for iDelay = 1:numel(vnDelay_fet)
    mr1_ = meanSubt_(trWav2(cvi1{iDelay},:,1));
    mr1_ = bsxfun(@rdivide, mr1_, sqrt(mean(mr1_.^2))); %zscore fast   
    tr1_ = repmat(mr1_, [1,1,nSites_spk]);
    tr2_ = meanSubt_(trWav2(cvi2{iDelay},:,:));
    cmrFet{iDelay} = permute(mean(tr1_ .* tr2_, 1), [3,2,1]);
end %for

if nargout==1
    varargout{1} = cell2mat_(cmrFet);
else
    varargout = cmrFet;
end
end %func


%--------------------------------------------------------------------------
function [trWav2, mrWav_ref] = spkwav_car_(trWav2, P, nSites_spk, viSite2_spk)
%function [trWav2, mrWav_ref] = trWav_car_sort_(trWav2, P)
%  trWav2: nT x nSpk x nSites, single
if nargin<3, nSites_spk = []; end
if nargin<4, viSite2_spk = []; end
vcSpkRef = get_set_(P, 'vcSpkRef', 'nmean');
if strcmpi(vcSpkRef, 'nmean')
    fSort_car = 1;
else
    fSort_car = -1; % no local subtraction
end
% fSort_car = get_set_(P, 'fSort_car', 1);
if isempty(nSites_spk)
    nSites_spk = 1 + 2 * P.maxSite - P.nSites_ref; % size(tnWav_spk, 2);
end
if nSites_spk==1, mrWav_ref=[]; return; end

% if P.nSites_ref==0, fSort_car = -1; end
% nSites_spk = 1 + 2 * P.maxSite; % size(tnWav_spk, 2);
dimm1 = size(trWav2);
[nT_spk, nSpk] = deal(dimm1(1), dimm1(2));
switch fSort_car
    case 1% use n sites having the least SD as reference sites
        if isempty(viSite2_spk)
            mrWav_ref = mean(trWav2(:,:,2:nSites_spk), 3);
        else
            trWav3 = trWav2(:,:,1:nSites_spk);
            trWav3(:,:,1) = 0;
            for iSpk1 = 1:numel(viSite2_spk)
                trWav3(:,iSpk1,viSite2_spk(iSpk1)) = 0;
            end
            mrWav_ref = sum(trWav3, 3) / (nSites_spk-2);
        end
    case 3
        mrWav_ref = mean(trWav2(:,:,2:end), 3);
    case 2
        [~, miSites_ref] = sort(squeeze_(var(trWav2)), 2, 'descend'); % use lest activities for ref
        %miSites_ref = miSites_ref(:,nSites_spk+1:end);
        miSites_ref = miSites_ref(:,3:nSites_spk);
        ti_dimm1 = repmat((1:nT_spk)', [1, nSpk, P.nSites_ref]);
        ti_dimm2 = repmat(1:nSpk, [nT_spk,1,P.nSites_ref]);
        ti_dimm3 = repmat(shiftdim(miSites_ref,-1), [nT_spk,1,1]);
        mrWav_ref = mean(trWav2(sub2ind(size(trWav2), ti_dimm1, ti_dimm2, ti_dimm3)),3);
    case 0
        mrWav_ref = mean(trWav2(:,:,nSites_spk+1:end), 3);
    case 4
        mrWav_ref = mean(trWav2(:,:,1:nSites_spk), 3);
    case 5
        mrWav_ref = mean(trWav2(:,:,2:end), 3);        
    case -1
        mrWav_ref = [];
end
trWav2 = trWav2(:,:,1:nSites_spk);
dimm2 = size(trWav2);
if ismatrix(trWav2), dimm2(end+1) = 1; end
if ~isempty(mrWav_ref)
    trWav2 = meanSubt_(reshape(bsxfun(@minus, reshape(trWav2,[],dimm2(3)), mrWav_ref(:)), dimm2));
%     trWav2 = (reshape(bsxfun(@minus, reshape(trWav2,[],dimm2(3)), mrWav_ref(:)), dimm2));
else
    trWav2 = meanSubt_(trWav2);
end
end %func


%--------------------------------------------------------------------------
function mr_ref = tr_sort_ref_(trWav_spk1, viSites_ref)
%  trWav_spk1: nT x nSpk x nSites, single
fUseMin = 0;
dimm1 = size(trWav_spk1);
tr0 = gather_(trWav_spk1);
%tr = permute(gather_(trWav_spk1), [3,1,2]);
mr_ref = zeros(dimm1([1,2]), 'like', tr0);
if fUseMin
    P = get0_('P');
    iT0 = 1 - P.spkLim(1);
    [~, miSites_ref] = sort(permute(tr0(iT0,:,:), [3,2,1]), 'ascend');
else % use std
    [~, miSites_ref] = sort(permute(var(tr0), [3,2,1]), 'descend'); % use lest activities for ref
end
miSites_ref = miSites_ref(viSites_ref,:);
for iSpk1=1:dimm1(2)
    %mr_ref(:,iSpk1) = mean(tr(miSites_ref(:,iSpk1),:,iSpk1));
    mr_ref(:,iSpk1) = mean(tr0(:,iSpk1,miSites_ref(:,iSpk1)), 3);
end
end %func


%--------------------------------------------------------------------------
function vr = mr2vr_sub2ind_(mr, vi1, vi2)
if isempty(mr), vr = []; return; end
if isempty(vi1), vi1 = 1:size(mr,1); end
if isempty(vi2), vi2 = 1:size(mr,2); end
vr = mr(sub2ind(size(mr), vi1(:), vi2(:)));
end %func


%--------------------------------------------------------------------------
function [fid, nBytes, header_offset] = fopen_(vcFile, vcMode)
if nargin < 2, vcMode = 'r'; end
try
    if matchFileExt_(vcFile, {'.ns5', '.ns2'})
        [fid, nBytes, header_offset] = fopen_nsx_(vcFile);
    else
        fid = fopen(vcFile, vcMode);
        nBytes = getBytes_(vcFile);   
        header_offset = 0;
    end
catch
    disperr_();
    fid = []; 
    nBytes = [];
end
end %func



%--------------------------------------------------------------------------
function nBytes = mem_max_(P)
try
    if P.fGpu
        S = gpuDevice(); % does not reset GPU
        nBytes = floor(S(1).AvailableMemory());
    else
        S = memory();
        nBytes = floor(S.MaxPossibleArrayBytes());
    end
catch
    nBytes = inf; % assume infinite memory
end
end %func


%--------------------------------------------------------------------------
function [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file, P)
% plan load file size according to the available memory and file size (nBytes_file1)
LOAD_FACTOR = 10; %GPU memory usage factor. 4x means 1/4 of GPU memory can be loaded

nSamples1 = floor(nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans);
% nSamples_max = floor(mem_max_(P) / P.nChans / 4); % Bound by MAX_BYTES_LOAD
if ~isfield(P, 'MAX_BYTES_LOAD'), P.MAX_BYTES_LOAD = []; end
if isempty(P.MAX_BYTES_LOAD), P.MAX_BYTES_LOAD = floor(mem_max_(P) / LOAD_FACTOR); end
if P.fft_thresh>0 %FFT memory limit    
    if isempty(P.MAX_LOAD_SEC)
        P.MAX_LOAD_SEC = 60; 
    else
        P.MAX_LOAD_SEC = min(P.MAX_LOAD_SEC, 60);
    end
    nSamples_max = min(floor(P.MAX_BYTES_LOAD / P.nChans), floor(P.sRateHz * P.MAX_LOAD_SEC)); % 1 min window
elseif isempty(P.MAX_LOAD_SEC)
   nSamples_max = floor(P.MAX_BYTES_LOAD / P.nChans / bytesPerSample_(P.vcDataType));
else
    nSamples_max = floor(P.sRateHz * P.MAX_LOAD_SEC);
end
nLoad1 = setlim_(ceil(nSamples1 / nSamples_max), [1, inf]); 
nSamples_load1 = min(nSamples1, nSamples_max);
if ~P.fTranspose_bin %load all in one, Catalin's format
    nLoad1 = 1; 
end
if nLoad1 == 1
    nSamples_load1 = nSamples1;
    nSamples_last1 = nSamples1;
else
    nSamples_last1 = mod(nSamples1, nSamples_load1);
    if nSamples_last1==0
        nSamples_last1 = nSamples_load1;
    elseif nSamples_last1 < nSamples_load1/2
        % if last part is too small increase the size
        nLoad1 = nLoad1 - 1;
        nSamples_last1 = nSamples_last1 + nSamples_load1;
    end
end
end %func


%--------------------------------------------------------------------------
% 17/9/13 JJJ: Created and tested
function vr = setlim_(vr, lim_)
% Set low and high limits
vr = min(max(vr, lim_(1)), lim_(2));
end %func


%--------------------------------------------------------------------------
function mnWav1 = load_file_preview_(fid_bin, P)
% preview
if P.nPad_filt > 0
    [mnWav1, vrWav_mean1, dimm_wav] = load_file_(fid_bin, P.nPad_filt, P);
    frewind_(fid_bin, dimm_wav, P.vcDataType);
else
    mnWav1 = [];
end
end %func


%--------------------------------------------------------------------------
function [mnWav1, vrWav_mean1, dimm_wav] = load_file_(fid_bin, nSamples_load1, P)
% load file to memory. output int16
% assume that the file is chan x time
% returns chans x time for efficient processing
% vrWav_mean1: average across chan output
if ischar(fid_bin)
    vcFile = fid_bin;
    [fid_bin, nBytes_file1] = fopen_(vcFile, 'r');
    nSamples_load1 = floor(nBytes_file1 / P.nChans / bytesPerSample_(P.vcDataType));    
else
    vcFile = [];
end
fSingle = 0; %output single
if P.fTranspose_bin
    dimm_wav = [P.nChans, nSamples_load1];    
else %Catalin's format
    dimm_wav = [nSamples_load1, P.nChans];
end
mnWav1 = fread_(fid_bin, dimm_wav, P.vcDataType);
[mnWav1, P.fGpu] = gpuArray_(mnWav1, P.fGpu);
switch(P.vcDataType)
    case 'uint16', mnWav1 = int16(single(mnWav1)-2^15);
    case {'single', 'double'}, mnWav1 = int16(mnWav1 / P.uV_per_bit);
end
if get_(P, 'fInverse_file'), mnWav1 = -mnWav1; end %flip the polarity

% extract channels
if P.fTranspose_bin
    if ~isempty(P.viSite2Chan), mnWav1 = mnWav1(P.viSite2Chan,:); end
    vrWav_mean1 = single(mean(mnWav1, 1)); %6x faster to transpose in dimm1
    if fSingle
        mnWav1 = single(mnWav1') * P.uV_per_bit; 
    else
        mnWav1 = mnWav1'; 
    end
else %Catalin's format. time x nChans
    if ~isempty(P.viSite2Chan), mnWav1 = mnWav1(:,P.viSite2Chan); end
    if ~isempty(P.tlim_load)
        nSamples = size(mnWav1,1);
        nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
        mnWav1 = mnWav1(nlim_load(1):nlim_load(end), :);
    end
    vrWav_mean1 = single(mean(mnWav1, 2)');
    if fSingle
        mnWav1 = single(mnWav1) * P.uV_per_bit; 
    end
end
if ~isempty(vcFile)
    fclose(fid_bin);
end
end %func


%--------------------------------------------------------------------------
function varargout = get_(varargin)
% retrieve a field. if not exist then return empty
% [val1, val2] = get_(S, field1, field2, ...)

if nargin==0, varargout{1} = []; return; end
S = varargin{1};
if isempty(S), varargout{1} = []; return; end

for i=2:nargin
    vcField = varargin{i};
    try
        varargout{i-1} = S.(vcField);
%         if ~isfield(S, vcField)
%             varargout{i-1} = [];
%         else
%             varargout{i-1} = S.(vcField);
%         end
    catch
        varargout{i-1} = [];
    end
end
end %func


%--------------------------------------------------------------------------
function [P, vcFile_prm] = loadParam_(vcFile_prm, fEditFile)

if nargin<2, fEditFile = 1; end
% Load prm file
if ~exist(vcFile_prm, 'file')
    error('.prm file does not exist: %s\n', vcFile_prm);
%     P=[]; return; 
end
P0 = file2struct_(jrcpath_(read_cfg_('default_prm', 0)));  %P = defaultParam();
P = file2struct_(vcFile_prm);
if ~isfield(P, 'template_file'), P.template_file = ''; end
if ~isempty(P.template_file)
    P0 = struct_merge_(P0, file2struct_(P.template_file));
end
P.vcFile_prm = vcFile_prm;
% todo: substitute bin file path
assert(isfield(P, 'vcFile'), sprintf('Check "%s" file syntax', vcFile_prm));

if ~exist(P.vcFile, 'file')
    P.vcFile = replacePath_(P.vcFile, vcFile_prm);
    if ~exist(P.vcFile, 'file')
        fprintf('vcFile not specified. Assuming multi-file format ''csFiles_merge''.\n');
    end
end
% Load prb file
if ~isfield(P, 'probe_file'), P.probe_file = P0.probe_file; end
try    
    probe_file_ = find_prb_(P.probe_file);
    if isempty(probe_file_)
        P.probe_file = replacePath_(P.probe_file, vcFile_prm); 
        if ~exist(P.probe_file, 'file'), error('prb file does not exist'); end
    else
        P.probe_file = probe_file_;
    end
    P0 = load_prb_(P.probe_file, P0);
catch
    fprintf('loadParam: %s not found.\n', P.probe_file);
end

% Load prb file
P = struct_merge_(P0, P);    

% check GPU
P.fGpu = ifeq_(license('test', 'Distrib_Computing_Toolbox'), P.fGpu, 0);
if P.fGpu, P.fGpu = ifeq_(gpuDeviceCount()>0, 1, 0); end

% Legacy support
if isfield(P, 'fTranspose'), P.fTranspose_bin = P.fTranspose; end

% Compute fields
P = struct_default_(P, 'fWav_raw_show', 0);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
P = struct_default_(P, 'vcFile_gt', '');
if isempty(P.vcFile_gt), P.vcFile_gt = subsFileExt_(P.vcFile_prm, '_gt.mat'); end
P.spkRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);
P.spkLim = round(P.spkLim_ms * P.sRateHz / 1000);
if isempty(P.spkLim_ms_fet), P.spkLim_ms_fet = P.spkLim_ms; end
% P.spkLim_fet = round(P.spkLim_ms_fet * P.sRateHz / 1000);
% P.slopeLim = round(P.slopeLim_ms * P.sRateHz / 1000);
if isempty(get_(P, 'nDiff_filt'))
    if isempty(get_(P, 'nDiff_ms_filt'))
        P.nDiff_filt = 0;
    else
        P.nDiff_filt = ceil(P.nDiff_ms_filt * P.sRateHz / 1000);
    end
end
if ~isfield(P, 'viShank_site'), P.viShank_site = []; end
try P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site); catch; end %find closest sites
P.sRateHz_lfp = P.sRateHz / P.nSkip_lfp;        %LFP sampling rate
P.bytesPerSample = bytesPerSample_(P.vcDataType);
P = struct_default_(P, 'vcFile_prm', subsFileExt_(P.vcFile, '.prm'));
if ~isempty(get_(P, 'gain_boost')), P.uV_per_bit = P.uV_per_bit / P.gain_boost; end
P.spkThresh = P.spkThresh_uV / P.uV_per_bit;
P = struct_default_(P, 'cvrDepth_drift', {});
P = struct_default_(P, {'maxSite_fet', 'maxSite_detect', 'maxSite_sort','maxSite_pix', 'maxSite_dip', 'maxSite_merge', 'maxSite_show'}, P.maxSite);
P = struct_default_(P, 'mrColor_proj', [.75 .75 .75; 0 0 0; 1 0 0]);
P.mrColor_proj = reshape(P.mrColor_proj(:), [], 3); %backward compatible
P = struct_default_(P, {'blank_thresh', 'thresh_corr_bad_site', 'tlim_load'}, []);
if isfield(P, 'rejectSpk_mean_thresh'), P.blank_thresh = P.rejectSpk_mean_thresh; end
P.vcFilter = get_filter_(P);
if isempty(get_(P, 'vcFilter_show'))
    P.vcFilter_show = P.vcFilter;
end
assert(validate_param_(P), 'Parameter file contains error.');
if fEditFile, edit(P.vcFile_prm); end % Show settings file
end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created and tested
function vcFilter = get_filter_(P)
if isfield(P, 'vcFilter')
    vcFilter = P.vcFilter;
else
    if get_(P, 'nDiff_filt') > 0
        vcFilter = 'sgdiff';
    else
        vcFilter = 'bandpass';
    end
end
end %func


%--------------------------------------------------------------------------
function [viSpk, vrSpk, viSite] = spikeMerge_(cviSpk, cvrSpk, P)
% provide spike index (cviSpk) and amplitudes (cvrSPk) per sites
% if nargin < 3
%     P = struct('spkLim', [-10 24], 'maxSite', 2.5);    
% end
% P.fParfor = 0; %debug

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
            disperr_();
        end
    end
    if ~P.fParfor
        fprintf('parfor disabled for spike merging.\n');
        for iSite = 1:nSites
            [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P);            
        end
    end
catch
    disperr_();
end

% merge parfor output and sort
viSpk = cell2mat_(cviSpkA);
vrSpk = cell2mat_(cvrSpkA);
viSite = cell2mat_(cviSiteA);
[viSpk, viSrt] = sort(viSpk); %sort by time
vrSpk = gather_(vrSpk(viSrt));
viSite = viSite(viSrt);
end %func


%--------------------------------------------------------------------------
function [viSpkA, vrSpkA, viSiteA] = spikeMerge_single_(viSpk, vrSpk, viSite, iSite, P)
% spkLim = int32(abs(P.spkLim(1))*[-1,1]); %user shorter lim
% fMerge_half = 1;
maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
nPad_pre = get_set_(P, 'nPad_pre', 0);
% refrac_factor = 2; %second refrac. set to 0 or [] to disable

nRefrac = int32(abs(P.spkRefrac));
spkLim = [-nRefrac, nRefrac];
% spkLim = int32(P.spkLim([1,end]));

% filter local sites only
vii1 = find(viSite == iSite);
viSpk1 = viSpk(vii1);
vrSpk1 = vrSpk(vii1);
viSiteNear = findNearSite_(P.mrSiteXY, iSite, maxDist_site_um);
vi2 = find(ismember(viSite, viSiteNear));
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
    if numel(vii2)==1, vlSpk1(iSpk1) = 1; continue; end %no other spikes detected
    [viSpk22, vrSpk22, viSite22] = select_vr_(viSpk2, vrSpk2, viSite2, vii2);

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

% apply spike merging on the same site
% if ~isempty(refrac_factor) && refrac_factor~=0
%     nRefrac2 = int32(double(nRefrac) * refrac_factor);
%     [viSpkA, vrSpkA, viSiteA] = spike_refrac_(viSpkA, vrSpkA, viSiteA, nRefrac2); %same site spikes
% end
end %func


%--------------------------------------------------------------------------
function viSiteNear = findNearSite_(mrSiteXY, iSite, maxDist_site_um)
vrDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
viSiteNear = find(vrDist <= maxDist_site_um);
end %func


%--------------------------------------------------------------------------
function [vlKeep_ref, vrMad_ref] = car_reject_(vrWav_mean1, P)
blank_period_ms = get_set_(P, 'blank_period_ms', 5); 
blank_thresh = get_set_(P, 'blank_thresh', []);
[vlKeep_ref, vrMad_ref] = deal([]);
% if isempty(blank_thresh) || blank_thresh==0, return; end

% tbin_ref = .01; %10 msec bin
vrWav_mean1 = single(vrWav_mean1);
nwin = round(P.sRateHz * blank_period_ms / 1000);
if nwin <= 1
    if nargout < 2
        vlKeep_ref = thresh_mad_(abs(vrWav_mean1), blank_thresh);
    else
        [vlKeep_ref, vrMad_ref] = thresh_mad_(abs(vrWav_mean1), blank_thresh);
    end
else
    vrRef_bin = std(reshape_vr2mr_(vrWav_mean1, nwin), 1,1);
    if nargout < 2
        vlKeep_ref = thresh_mad_(vrRef_bin, blank_thresh);
    else
        [vlKeep_ref, vrMad_ref] = thresh_mad_(vrRef_bin, blank_thresh);        
        vrMad_ref = expand_vr_(vrMad_ref, nwin, size(vrWav_mean1));
    end
    vlKeep_ref = expand_vr_(vlKeep_ref, nwin, size(vrWav_mean1));
end
% figure; plot(vrMad_ref); hold on; plot(find(~vlKeep_ref), vrMad_ref(~vlKeep_ref), 'r.')
end %func


%--------------------------------------------------------------------------
% 8/16/17 JJJ: created and tested
function vr1 = expand_vr_(vr, nwin, dimm1);
if nargin<3, dimm1 = [numel(vr) * nwin, 1]; end
if islogical(vr)
    vr1 = false(dimm1);
else
    vr1 = zeros(dimm1, 'like', vr);
end
vr = repmat(vr(:)', [nwin, 1]);
vr = vr(:);
[n,n1] = deal(numel(vr), numel(vr1));
if n1 > n
    vr1(1:n) = vr;
    vr1(n+1:end) = vr1(n);
elseif numel(vr1) < n
    vr1 = vr(1:n1);
else
    vr1 = vr;
end
end %func


%--------------------------------------------------------------------------
function flag = is_detected_(P)
% return true if already detected. .spkwav file must exist
flag = exist(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), 'file');
end %func


%--------------------------------------------------------------------------
function flag = is_sorted_(P)
% return true if already detected. .spkwav file must exist
S0 = load0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
flag = ~isempty(get0_('S_clu'));
% flag = exist(strrep(P.vcFile_prm, '.prm', '_clu.mat'), 'file');
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Added prb directory
function P = load_prb_(vcFile_prb, P)
% append probe file to P

% Find the probe file
vcFile_prb = find_prb_(vcFile_prb);
if isempty(vcFile_prb)
    error(['Probe file does not exist: ', vcFile_prb]);
end

P.probe_file = vcFile_prb;
%     [P.viSite2Chan, P.mrSiteXY, P.vrSiteHW, P.cviShank] = read_prb_file(vcFile_prb);
S_prb = file2struct_(vcFile_prb);
P.viSite2Chan = S_prb.channels;
P.mrSiteXY = S_prb.geometry;
P.vrSiteHW = S_prb.pad;
if ~isfield(S_prb, 'shank')
    P.viShank_site = ones(size(S_prb.channels)); 
else
    P.viShank_site = S_prb.shank;
end
S_prb = remove_struct_(S_prb, 'channels', 'geometry', 'pad', 'ref_sites', ...
    'viHalf', 'i', 'vcFile_file2struct', 'shank');

% P = copyStruct_(P, S_prb, {'cviShank', 'maxSite', 'um_per_pix'});
P.viChan_aux = setdiff(1:P.nChans, 1:max(P.viSite2Chan)); %aux channel. change for
P = struct_merge_(P, S_prb);
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile = search_file_(vcFile, csDir)
% Search file in the provided directory location if it doesn't exist
if ischar(csDir), csDir = {csDir}; end

if ~exist(vcFile, 'file')
    for iDir = 1:numel(csDir)
        vcFile_ = subsDir_(vcFile, csDir{iDir});
        if exist(vcFile_, 'file')
            vcFile = vcFile_;
            return;
        end
    end
    vcFile = []; % file not found
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function vcFile_new = subsDir_(vcFile, vcDir_new)
% Substitute dir
[vcDir, vcFile, vcExt] = fileparts(vcFile);
vcFile_new = fullfile(vcDir_new, [vcFile, vcExt]);
end % func


%--------------------------------------------------------------------------
function S = remove_struct_(S, varargin)
% remove fields from a struct
for i=1:numel(varargin)
    if isfield(S, varargin{i})
        S = rmfield(S, varargin{i});
    end
end
end %func


%--------------------------------------------------------------------------
% 9/13/17 JJJ: Edge case error fixed (vi2 indexing) for small block many loads
% 8/17/17 JJJ: fix for missing spikes
function viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min)
% nneigh_min: number of neighbors around the spike below the threshold
%  0,1,2. # neighbors of minimum point below negative threshold 
% thresh1: absolute value. searching for negative peaks only

if nargin<3, nneigh_min = []; end
if isempty(nneigh_min), nneigh_min = 1; end

viSpk1 = [];
if isempty(vrWav1), return; end
vl1 = vrWav1 < -abs(thresh1);
vi2 = find(vl1);
%vi2 = find(vrWav1 < -thresh1);
if isempty(vi2), return; end

if vi2(1)<=1
    if numel(vi2) == 1, return; end
    vi2(1) = []; 
end    
if vi2(end)>=numel(vrWav1)
    if numel(vi2) == 1, return; end
    vi2(end) = []; 
end
vrWav12 = vrWav1(vi2);
viSpk1 = vi2(vrWav12 <= vrWav1(vi2+1) & vrWav12 <= vrWav1(vi2-1));
if isempty(viSpk1), return; end
% viSpk1 = vi2(find(diff(diff(vrWav1(vi2))>0)>0) + 1); % only negative peak

% if viSpk1(1) <= 1, viSpk1(1) = 2; end
% if viSpk1(end) >= numel(vrWav1), viSpk1(end) = numel(vrWav1)-1; end
switch nneigh_min
    case 1
        viSpk1 = viSpk1(vl1(viSpk1-1) | vl1(viSpk1+1));
%         viSpk1 = viSpk1(vrWav1(viSpk1-1) < -thresh1 | vrWav1(viSpk1+1) < -thresh1);
    case 2
        viSpk1 = viSpk1(vl1(viSpk1-1) & vl1(viSpk1+1));
%         viSpk1 = viSpk1(vrWav1(viSpk1-1) < -thresh1 & vrWav1(viSpk1+1) < -thresh1);
end
end %func


%--------------------------------------------------------------------------
function vr = cell2mat_(cvr)
% remove empty
vi = find(cellfun(@(x)~isempty(x), cvr));
vr = cell2mat(cvr(vi));
end %func


%--------------------------------------------------------------------------
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
    if ia==0
        if vi(i) >= a, ia = i; end
    end
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


%--------------------------------------------------------------------------
function dimm_mr = write_bin_(vcFile, mr)
t1=tic;
dimm_mr = size(mr);
if isempty(mr), return; end
if isstruct(mr)
    save(vcFile, '-struct', 'mr', '-v7.3');
else
    fidw = fopen(vcFile, 'W');
    fwrite(fidw, mr, class(mr));
    fclose(fidw);
end
fprintf('Writing to %s took %0.1fs\n', vcFile, toc(t1));
end %func


%--------------------------------------------------------------------------
function write_struct_(vcFile, S)
% Write a struct S to file vcFile
try
    warning off
    t1=tic;    
%     S = struct_remove_handles(S); %remove figure handle
    save(vcFile, '-struct', 'S');
    fprintf('Wrote to %s, took %0.1fs.\n', vcFile, toc(t1));
catch
    fprintf(2, 'Writing struct to file %s failed.\n', vcFile);
end
end %func


%--------------------------------------------------------------------------
function mnWav = load_bin_(vcFile, vcDataType, dimm)
% mnWav = load_bin_(vcFile, dimm, vcDataType)
% mnWav = load_bin_(fid, dimm, vcDataType)

if nargin<2, vcDataType = []; end
if nargin<3, dimm = []; end
if isempty(vcDataType), vcDataType = 'int16'; end
if ischar(vcFile)
    fid = []; 
    if ~exist(vcFile, 'file')
        fprintf(2, 'File does not exist: %s\n', vcFile);
        mnWav=[]; 
        return; 
    end
    if isempty(dimm)
        S_file = dir(vcFile);
        nData = S_file(1).bytes / bytesPerSample_(vcDataType);
        dimm = [nData, 1]; %return column
    end
else % fid directly passed
    fid = vcFile;
    vcFile = [];
    if isempty(dimm), dimm = inf; end
end
try
    t1 = tic;
    if isempty(fid), fid = fopen(vcFile, 'r'); end
    mnWav = fread_(fid, dimm, vcDataType);
    if ~isempty(vcFile), fclose(fid); end %close fid if file name passed
    if ~isempty(vcFile), fprintf('Loading %s took %0.1fs\n', vcFile, toc(t1)); end
catch
    disperr_();
    mnWav = [];
end
end %func


%--------------------------------------------------------------------------
function mnWav = load_bin_multi_(fid, cvi_lim_bin, P)
mnWav = cell(size(cvi_lim_bin));
fpos = 0;
for i=1:numel(cvi_lim_bin)
    lim1 = cvi_lim_bin{i};
    if i>1, fseek_(fid, lim1(1), P); end    
    mnWav{i} = load_bin_(fid, P.vcDataType, [P.nChans, diff(lim1)+1]);
    if i==1, fpos = ftell(fid); end
end %for
mnWav = cell2mat_(mnWav);
fseek(fid, fpos, 'bof'); % restore the file position
end %func


%--------------------------------------------------------------------------
function [vrY, vrY2] = centroid_mr_(mrVpp, vrYe, mode1)
% [vrX, vrY] = centroid_mr_(mrVpp, vrPos)
% [mrXY] = centroid_mr_(mrVpp, vrPos)
% mrVpp: nSites x nSpk

if nargin<3, mode1 = 1; end
if isrow(vrYe), vrYe = vrYe'; end
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


%--------------------------------------------------------------------------
function S_clu = fet2clu_(S0, P)
% can process different shanks separately
global trFet_spk
fRepeat_clu = get_set_(P, 'fRepeat_clu', 1);
fprintf('Clustering\n');

S_clu = cluster_spacetime_(S0, P);        
S_clu = postCluster_(S_clu, P);

if fRepeat_clu
    trFet_spk0 = trFet_spk;
    nSites_fet = P.maxSite*2+1-P.nSites_ref;
    nFetPerSite = size(trFet_spk,1) / nSites_fet;
    vrSnr_clu = S_clu_snr_(S_clu);
    %         vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 2^(-iRepeat_clu+1)); %exponential selection
    vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/2); %ilnear selection
    vlRedo_spk = ismember(S_clu.viClu, find(vlRedo_clu));        

    % top half
%     S_clu_A = postCluster_(cluster_spacetime_(S0, P, ~vlRedo_spk), P);

    % Bottom half (reduce feature diemnsion)    
    mlFet_ = false(nSites_fet, nFetPerSite);
    nSites_fet = ceil(nSites_fet*.75);
    mlFet_(1:nSites_fet, :) = 1;
    trFet_spk = trFet_spk(find(mlFet_),:,:);
    S_clu_B = postCluster_(cluster_spacetime_(S0, P, vlRedo_spk), P);        

    % combine
%     S_clu.viClu(~vlRedo_spk) = S_clu_A.viClu;
%     S_clu.viClu(vlRedo_spk) = S_clu_B.viClu + max(S_clu_A.viClu);
    S_clu = S_clu_combine_(S_clu, S_clu_B, vlRedo_clu, vlRedo_spk);
    trFet_spk = trFet_spk0; %restore
end
% S_clu = clu2wav_(S_clu, S0.viSite_spk, tnWav_spk, tnWav_raw);
S_clu = post_merge_(S_clu, P, 0);
S_clu.viClu_auto = S_clu.viClu;
fprintf('\tClustering took %0.1f s\n', S_clu.t_runtime);
end %func


%--------------------------------------------------------------------------
function [vi1, vl1] = partition_vi_(vi, n, nBins, iBin)
lim1 = round([iBin-1, iBin] / nBins * n) + 1;
lim1 = min(max(lim1, 1), n+1);
vl1 = vi >= lim1(1) & vi < lim1(2);
vi1 = vi(vl1);
end %func


%--------------------------------------------------------------------------
function d = eucl2_dist_(X, Y)
% a: m x d1; b: m x d2
% aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b; 
% d = sqrt(abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab));
% X = [mrFet1_, mrFet2_];
d = bsxfun(@plus, sum(Y.^2), bsxfun(@minus, sum(X.^2)', 2*X'*Y));
end %func


%--------------------------------------------------------------------------
function S_clu = cluster_spacetime_(S0, P, vlRedo_spk)
% this clustering is natively suited for 2D electrode arrays and drifting dataset
% There is no x,y position in the clustering dataset
global trFet_spk
if ~isfield(P, 'CHUNK'), P.CHUNK = 16; end
if ~isfield(P, 'fTwoStep'), P.fTwoStep = 0; end
if ~isfield(P, 'mrSiteXY'), P.mrSiteXY = []; end
if ~isfield(P, 'min_count'), P.min_count = []; end
if nargin<3, vlRedo_spk=[]; end

% g = gpuDevice();
t_func = tic;
nSites = numel(P.viSite2Chan);
nSpk = numel(S0.viTime_spk);
vrRho = zeros(nSpk, 1, 'single');
vrDelta = zeros(nSpk, 1, 'single');   
viNneigh = zeros(nSpk, 1, 'uint32');
vrDc2_site = zeros(nSites, 1, 'single');
nTime_clu = get_set_(P, 'nTime_clu', 1);
P.nTime_clu = nTime_clu;
P.dc_subsample = 1000; 
P.fGpu = 1;

% clear memory
cuda_rho_();
cuda_delta_();
if get_set_(P, 'fDc_global', 0)
    dc2 = calc_dc2_(S0, P, vlRedo_spk);
else
    dc2 = [];
end

%-----
% Calculate Rho
fprintf('Calculating Rho\n\t'); t1=tic;
for iSite = 1:nSites
    [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk);
    if isempty(mrFet12_), continue; end
    if isempty(dc2)
        dc2_ = compute_dc2_(mrFet12_, viiSpk12_ord_, n1_, n2_, P); % Compute DC in CPU
    else
        dc2_ = dc2.^2;
    end
    [mrFet12_, viiSpk12_ord_] = multifun_(@gpuArray_, mrFet12_, viiSpk12_ord_);
    vrRho_ = cuda_rho_(mrFet12_, viiSpk12_ord_, n1_, n2_, dc2_, P); 
    viSpk_site_ = S0.cviSpk_site{iSite};
    if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
    vrRho(viSpk_site_) = gather_(vrRho_);
    vrDc2_site(iSite) = dc2_;
    [mrFet12_, viiSpk12_ord_, vrRho_] = deal([]);
    fprintf('.');    
end
% [vrRho, vrDc2_site] = gather_(vrRho, vrDc2_site);
fprintf('\n\ttook %0.1fs\n', toc(t1));

%-----
% Calculate Delta
fprintf('Calculating Delta\n\t'); t2=tic;
for iSite = 1:nSites
    [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk);
    if isempty(mrFet12_), continue; end    
    viiRho12_ord_ = rankorder_(vrRho(viSpk12_), 'descend');    
    [mrFet12_, viiRho12_ord_, viiSpk12_ord_] = ...
        multifun_(@gpuArray_, mrFet12_, viiRho12_ord_, viiSpk12_ord_);
    try
        [vrDelta_, viNneigh_] = cuda_delta_(mrFet12_, viiSpk12_ord_, viiRho12_ord_, n1_, n2_, vrDc2_site(iSite), P);
        [vrDelta_, viNneigh_] = gather_(vrDelta_, viNneigh_);
    catch
        disperr_(sprintf('error at site# %d', iSite));
    end
    viSpk_site_ = S0.cviSpk_site{iSite};
    if ~isempty(vlRedo_spk), viSpk_site_ = viSpk_site_(vlRedo_spk(viSpk_site_)); end
    vrDelta(viSpk_site_) = vrDelta_;
    viNneigh(viSpk_site_) = viSpk12_(viNneigh_);
    [mrFet12_, viiRho12_ord_, viiSpk12_ord_] = deal([]); %vrDelta_, viNneigh_, viSpk12_
    fprintf('.');
end
% Deal with nan delta
viNan_delta = find(isnan(vrDelta));
if ~isempty(viNan_delta)
    vrDelta(viNan_delta) = max(vrDelta);
end
% [vrDelta, viNneigh] = multifun_(@gather_, vrDelta, viNneigh);
fprintf('\n\ttook %0.1fs\n', toc(t2));

if ~isempty(vlRedo_spk)
    vrRho = vrRho(vlRedo_spk);
    vrDelta = vrDelta(vlRedo_spk);        
    viNneigh = reverse_lookup_(viNneigh(vlRedo_spk), find(vlRedo_spk));
end

%-----
% package
% if P.fGpu
%     [vrRho, vrDelta, viNneigh] = multifun_(@gather_, vrRho, vrDelta, viNneigh);
% end
t_runtime = toc(t_func);
trFet_dim = size(trFet_spk); %[1, size(mrFet1,1), size(mrFet1,2)]; %for postCluster
[~, ordrho] = sort(vrRho, 'descend');
S_clu = struct('rho', vrRho, 'delta', vrDelta, 'ordrho', ordrho, 'nneigh', viNneigh, ...
    'P', P, 't_runtime', t_runtime, 'halo', [], 'viiSpk', [], 'trFet_dim', trFet_dim, 'vrDc2_site', vrDc2_site);

% figure; loglog(vrRho, vrDelta, '.');
end %func


%--------------------------------------------------------------------------
function vrRho1 = cuda_rho_(mrFet12, viiSpk12_ord, n1, n2, dc2, P)
% Ultimately use CUDA to do this distance computation
% mrFet12_: already in GPU
% viiSpk12: ordered list
persistent CK nC_
if nargin==0, nC_ = 0; return; end
if isempty(nC_), nC_ = 0; end
[nC, n12] = size(mrFet12); %nc is constant with the loop
dn_max = int32(round((n1+n2) / P.nTime_clu));
nC_max = get_set_(P, 'nC_max', 45);
dc2 = single(dc2);
FLAG_FIXN = 0; %flag variable number of neighbors (otherwise fixed to 2*dn_max+1)
if P.fGpu && FLAG_FIXN == 0
    try        
        if (nC_ ~= nC) % create cuda kernel
            nC_ = nC;
            CK = parallel.gpu.CUDAKernel('jrc3_cuda_rho.ptx','jrc3_cuda_rho.cu');
            CK.ThreadBlockSize = [P.nThreads, 1];          
            CK.SharedMemorySize = 4 * P.CHUNK * (2 + nC_max + 2 * P.nThreads); % @TODO: update the size
        end
        CK.GridSize = [ceil(n1 / P.CHUNK / P.CHUNK), P.CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
        vrRho1 = zeros([1, n1], 'single', 'gpuArray'); 
        vnConst = int32([n1, n12, nC, dn_max, get_set_(P, 'fDc_spk', 0)]);
        vrRho1 = feval(CK, vrRho1, mrFet12, viiSpk12_ord, vnConst, dc2);
        return;
    catch        
        disperr_();
        nC_ = 0;
        fprintf(2, 'CUDA kernel failed. Re-trying in CPU\n');
    end
end

viiSpk1_ord= viiSpk12_ord(1:n1)';
mlKeep12_ = abs(bsxfun(@minus, viiSpk12_ord, viiSpk1_ord)) <= dn_max;
if FLAG_FIXN == 0
    vrRho1 = sum(eucl2_dist_(mrFet12, mrFet12(:,1:n1)) < dc2 & mlKeep12_) - 1; %do not include self    
    vrRho1 = single(vrRho1 ./ sum(mlKeep12_));
else
    mlKeep12_(1:min(2*dn_max+1,n12), viiSpk1_ord <= min(dn_max,n1)) = 1;
    mlKeep12_(max(n12-2*dn_max,1):n12,  viiSpk1_ord >= max(n1-dn_max,1)) = 1;
    vrRho1 = sum(eucl2_dist_(mrFet12, mrFet12(:,1:n1)) < dc2 & mlKeep12_) - 1; %do not include self    
    vrRho1 = single(vrRho1 ./ (2*dn_max));
end
end


%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1] = cuda_delta_(mrFet12, viiSpk12_ord, viiRho12_ord, n1, n2, dc2, P)
% Ultimately use CUDA to do this distance computation
persistent CK nC_
if nargin==0, nC_ = 0; return; end
if isempty(nC_), nC_ = 0; end
[nC, n12] = size(mrFet12); %nc is constant with the loop
dn_max = int32(round((n1+n2) / P.nTime_clu));
nC_max = get_set_(P, 'nC_max', 45);
if P.fGpu
    try
        if (nC_ ~= nC) % create cuda kernel
            nC_ = nC;
            CK = parallel.gpu.CUDAKernel('jrc3_cuda_delta.ptx','jrc3_cuda_delta.cu');
            CK.ThreadBlockSize = [P.nThreads, 1];          
            CK.SharedMemorySize = 4 * P.CHUNK * (3 + nC_max + 2*P.nThreads); % @TODO: update the size
        end
        CK.GridSize = [ceil(n1 / P.CHUNK / P.CHUNK), P.CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]    
        vrDelta1 = zeros([1, n1], 'single', 'gpuArray'); 
        viNneigh1 = zeros([1, n1], 'uint32', 'gpuArray'); 
        vnConst = int32([n1, n12, nC, dn_max, get_set_(P, 'fDc_spk', 0)]);
        [vrDelta1, viNneigh1] = feval(CK, vrDelta1, viNneigh1, mrFet12, viiSpk12_ord, viiRho12_ord, vnConst, dc2);
        % [vrDelta1_, viNneigh1_] = deal(vrDelta1, viNneigh1);
        return;
    catch        
        disperr_();
        nC_ = 0;
        fprintf(2, 'CUDA kernel failed. Re-trying in CPU\n');
    end
end

mrDist12_ = eucl2_dist_(mrFet12, mrFet12(:,1:n1));  %not sqrt
mlRemove12_ = bsxfun(@ge, viiRho12_ord, viiRho12_ord(1:n1)') ...
    | abs(bsxfun(@minus, viiSpk12_ord, viiSpk12_ord(1:n1)')) > dn_max;
mrDist12_(mlRemove12_) = nan;
[vrDelta1, viNneigh1] = min(mrDist12_);
end
    
    
%--------------------------------------------------------------------------
function dc2_ = compute_dc2_(mrFet12, viiSpk12_ord, n1_, n2_, P)
% subsample
fSubsample12 = 1;
if get_set_(P, 'fDc_spk', 0)
    dc2_ = (P.dc_percent/100).^2; 
    return; 
end % spike-specific dc

vi1_ = subsample_vr_(1:n1_, P.dc_subsample);
viiSpk1_ord_ = viiSpk12_ord(vi1_);
mrFet1_ = mrFet12(:,vi1_);

if fSubsample12
%     vi12_ = subsample_vr_(1:(n1_+n2_), P.dc_subsample * 4);
    vi12_ = subsample_vr_(1:n1_, P.dc_subsample * 4);
    viiSpk12_ord = viiSpk12_ord(vi12_);
    mrFet12 = mrFet12(:,vi12_);
end

mrDist11_2_ = eucl2_dist_(mrFet12, mrFet1_);
mlKeep11_2_ = abs(bsxfun(@minus, viiSpk12_ord, viiSpk1_ord_')) < (n1_+n2_) / P.nTime_clu;
mrDist11_2_(~mlKeep11_2_) = nan;
mrDist_sub = subsample_mr_(mrDist11_2_, P.dc_subsample, 2);
mrDist_sub(mrDist_sub<=0) = nan;
dc2_ = nanmedian(quantile(mrDist_sub, P.dc_percent/100));  
if isnan(dc2_)
    dc2_ = quantile(mrDist_sub(:), P.dc_percent/100);
end
end %func


%--------------------------------------------------------------------------
function [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk)
% decide whether to use 1, 2, or 3 features
if nargin<5, vlRedo_spk = []; end
nFet_use = get_set_(P, 'nFet_use', 2);
[viSpk1_, viSpk2_] = multifun_(@int32, S0.cviSpk_site{iSite}, S0.cviSpk2_site{iSite});
if ~isempty(vlRedo_spk)
    viSpk1_ = viSpk1_(vlRedo_spk(viSpk1_));
    viSpk2_ = viSpk2_(vlRedo_spk(viSpk2_));
end
if isempty(viSpk1_)
    [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = deal([]);
    return;
end
if isempty(viSpk2_), nFet_use=1; end
switch nFet_use
    case 3
        viSpk3_ = int32(S0.cviSpk3_site{iSite});
        if ~isempty(viSpk3_), viSpk3_ = viSpk3_(vlRedo_spk(viSpk3_)); end
        mrFet12_ = [squeeze_(trFet_spk(:,1,viSpk1_)), squeeze_(trFet_spk(:,2,viSpk2_)), squeeze_(trFet_spk(:,3,viSpk3_))];
        viSpk12_ = [viSpk1_; viSpk2_; viSpk3_];
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_) + numel(viSpk3_));
    case 2
        mrFet12_ = [squeeze_(trFet_spk(:,1,viSpk1_)), squeeze_(trFet_spk(:,2,viSpk2_))];
        viSpk12_ = [viSpk1_; viSpk2_];
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
    case 1
        mrFet12_ = squeeze_(trFet_spk(:,1,viSpk1_));
        viSpk12_ = viSpk1_;
        [n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
end
if get_set_(P, 'fSqrt_fet', 0), mrFet12_ = signsqrt_(mrFet12_); end
if get_set_(P, 'fLog_fet', 0), mrFet12_ = signlog_(mrFet12_); end
viiSpk12_ord_ = rankorder_(viSpk12_, 'ascend');
end %func


%--------------------------------------------------------------------------
function [mrFet1_, mrFet2_, viSpk1_, viSpk2_, n1_, n2_, viiSpk12_ord_] = fet12_site_1_(mrFet1, mrFet2, S0, P, iSite)

[viSpk1_, viSpk2_] = deal(S0.cviSpk_site{iSite}, S0.cviSpk2_site{iSite});
[mrFet1_, mrFet2_] = deal(mrFet1(:,viSpk1_), mrFet2(:,viSpk2_));
[n1_, n2_] = deal(numel(viSpk1_), numel(viSpk2_));
if P.fGpu
   [mrFet1_, mrFet2_, viSpk1_, viSpk2_] = multifun_(@gpuArray_, mrFet1_, mrFet2_, viSpk1_, viSpk2_);
end
viiSpk12_ord_ = rankorder_([viSpk1_; viSpk2_], 'ascend');
end %func


%--------------------------------------------------------------------------
function validate_(P)
global tnWav_spk tnWav_gt S_gt vcFile_prm_
% S0 = load_cached_(P, 0);
fMergeCheck = 0; %kilosort-style validation
snr_thresh_score = 10;
snr_thresh_stat = 6;

S0 = load_(strrep(P.vcFile_prm, '.prm', '_jrc.mat')); 
if isempty(tnWav_spk)
    tnWav_spk = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), P.vcDataType, S0.dimm_spk);
end
S_clu = S0.S_clu;

% Load ground truth file
if ~exist(P.vcFile_gt, 'file'), P.vcFile_gt = subsFileExt_(P.vcFile, '_gt.mat'); end
if ~strcmpi(vcFile_prm_, P.vcFile_prm), vcFile_prm_=[]; end
if isempty(vcFile_prm_) % cache ground truth result
    vcFile_prm_ = P.vcFile_prm;
    S_gt = load_gt_(P.vcFile_gt, P);
    if isempty(S_gt), fprintf(2, 'Groundtruth does not exist. Run "jrclust import" to create a groundtruth file.\n'); return; end
    [S_gt, tnWav_gt] = gt2spk_(S_gt, P, snr_thresh_stat);
end
S_score = struct(...
    'vrVmin_gt', S_gt.vrVmin_clu, 'vnSite_gt', S_gt.vnSite_clu, ...
    'vrSnr_gt', S_gt.vrSnr_clu, 'vrSnr_min_gt', S_gt.vrSnr_clu, ...
    'trWav_gt', S_gt.trWav_clu, 'viSite_gt', S_gt.viSite_clu);
S_score.cviSpk_gt = S_gt.cviSpk_clu;

% Compare S_clu with S_gt
nSamples_jitter = round(P.sRateHz / 1000); %1 ms jitter
fprintf('verifying cluster...\n'); 
[mrMiss, mrFp, vnCluGt, miCluMatch, S_score_clu] = ...
    clusterVerify(S_gt.viClu, S_gt.viTime, S_clu.viClu, S0.viTime_spk, nSamples_jitter);  %S_gt.viTime
% viClu_spk = S_score_clu.cviHit_gt
if fMergeCheck, compareClustering2_(S_gt.viClu, S_gt.viTime, S_clu.viClu+1, S0.viTime_spk); end

Sgt = S_gt; %backward compatibility
S_score = struct_add_(S_score, mrMiss, mrFp, vnCluGt, miCluMatch, P, Sgt, S_score_clu);
S_score.cviTime_clu = S_clu.cviSpk_clu(S_score_clu.viCluMatch)';
S_score.vrVrms_site = single(S0.vrThresh_site) / S0.P.qqFactor;

fprintf('SNR_gt (Vp/Vrms): %s\n', sprintf('%0.1f ', S_score.vrSnr_gt));
fprintf('nSites>thresh (GT): %s\n', sprintf('%d ', S_score.vnSite_gt));

write_struct_(strrep(P.vcFile_prm, '.prm', '_score.mat'), S_score);

set0_(S_score);       
assignWorkspace_(S_score); %put in workspace

figure; set(gcf,'Name',P.vcFile_prm);
subplot 121; plot_cdf_(S_score.S_score_clu.vrFp); hold on; plot_cdf_(S_score.S_score_clu.vrMiss); 
legend({'False Positives', 'False Negatives'}); ylabel('CDF'); grid on; xlabel('Cluster count');

subplot 122; hold on;
plot(S_score.vrSnr_min_gt, S_score.S_score_clu.vrFp, 'b.', S_score.vrSnr_min_gt, S_score.S_score_clu.vrMiss, 'r.');
legend({'False Positives', 'False Negatives'}); ylabel('score'); grid on; xlabel('SNR (Vp/Vrms)');

disp_score_(S_score.vrSnr_min_gt, S_score.S_score_clu.vrFp, S_score.S_score_clu.vrMiss, snr_thresh_score);

end %func


%--------------------------------------------------------------------------
function S_gt = load_gt_(vcFile_gt, P)
% S_gt contains viTime and viClu
if nargin<2, P = get0_('P'); end
if ~exist(vcFile_gt, 'file'), S_gt=[]; return; end
S = load(vcFile_gt);
if isfield(S, 'S_gt')
    S_gt = S.S_gt;
elseif isfield(S, 'viClu') && isfield(S, 'viTime')
    S_gt = S;
elseif isfield(S, 'viClu') && isfield(S, 'viSpk')
    S_gt.viTime = S.viSpk;    
    S_gt.viClu = S.viClu;    
else
    % Convert Nick's format to JRCLUST fomat
    if isfield(S, 'gtTimes')
        S_gt.viTime = cell2mat_(S.gtTimes');
        S_gt.viClu = cell2mat_(arrayfun(@(i)ones(size(S.gtTimes{i}))*i, 1:numel(S.gtTimes), 'UniformOutput', 0)');
    else
        error('no field found.');
    end
    [S_gt.viTime, ix] = sort(S_gt.viTime, 'ascend');
    S_gt.viClu = S_gt.viClu(ix);
end
if ~isempty(get_(P, 'tlim_load'))
    nSamples = double(S_gt.viTime(end));
    nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
    viKeep = find(S_gt.viTime >= nlim_load(1) & S_gt.viTime <= nlim_load(2));
    [S_gt.viTime, S_gt.viClu] = multifun_(@(x)x(viKeep), S_gt.viTime, S_gt.viClu);
end
[viClu_unique, ~, viClu] = unique(S_gt.viClu);
if max(S_gt.viClu) > numel(viClu_unique)
    S_gt.viClu = viClu;
end
end %func


%--------------------------------------------------------------------------
function [S_gt, tnWav_spk, tnWav_raw] = gt2spk_(S_gt, P, snr_thresh)
% convert ground truth to spike waveforms
fSubtract_nmean = 0;
% fSubtract_ref = 1;
% P.fGpu = 0;
MAX_SAMPLE = 1000; %mean calc
nSubsample_clu = get_set_(P, 'nSubsample_clu', 8);

if nargin<3, snr_thresh = []; end

fProcessRaw = (nargout == 3);
t1 = tic;
fprintf('Computing ground truth units...\n');
viClu = int32(S_gt.viClu);
viTime_spk = int32(S_gt.viTime);
nSites = numel(P.viSite2Chan);

% load entire raw waveform to memory
if ~P.fTranspose_bin
    [mnWav, vrWav_mean] = load_file_(P.vcFile, [], P);
    mnWav = int16(fft_clean_(mnWav, P));
    if fProcessRaw
        tnWav_raw = permute(mn2tn_gpu_(mnWav, P.spkLim * 2, viTime_spk), [1,3,2]);
    end
    [mnWav, ~] = filt_car_(mnWav, P);
%     mnWav = mnWav_filt_(mnWav, P); % Apply filtering in RAM
    if fSubtract_nmean % Apply nmean CAR to ground truth spikes (previous standard)
        P1=P; P1.vcCommonRef = 'nmean'; mnWav = wav_car_(mnWav, P1);
    end
    tnWav_spk = permute(mn2tn_gpu_(mnWav, P.spkLim, viTime_spk), [1,3,2]);
    vrVrms_site = gather_(mr2rms_(mnWav, 1e5));
    clear mnWav;
else % real ground truth: must block load and filter.    
    viClu = viClu(1:nSubsample_clu:end);
    viTime_spk = viTime_spk(1:nSubsample_clu:end);
    [tnWav_spk, vrVrms_site] = file2spk_gt_(P, viTime_spk);
    snr_thresh = []; %no SNR threshold
end

% determine mean spikes
nClu = max(viClu);
trWav_clu = zeros(size(tnWav_spk,1), nSites, nClu, 'single');
if fProcessRaw
    trWav_raw_clu = zeros(size(tnWav_raw,1), nSites, nClu, 'single');
else
    trWav_raw_clu = [];
end
cviSpk_clu = arrayfun(@(iClu)int32(find(viClu == iClu)), 1:nClu, 'UniformOutput', 0);
for iClu=1:nClu
    viSpk_clu1 = cviSpk_clu{iClu};
%     viSpk_clu1 = viSpk_clu1(
    viSpk1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
    if isempty(viSpk1), continue; end
    trWav_clu(:,:,iClu) = mean(tnWav_spk(:,:,viSpk1), 3); %multiply by scaling factor?
    if fProcessRaw
        trWav_raw_clu(:,:,iClu) = mean_tnWav_raw_(tnWav_raw(:,:,viSpk1), P);
    end
    fprintf('.');
end

% if fSubtract_ref
% % apply surrounding reference subtraction
%     trWav_clu = spkwav_car_(trWav_clu, spkwav_car_init_(P));
% end

% Find center location and spike SNR
mrVmin_clu = shiftdim(min(trWav_clu));
[vrVmin_clu, viSite_clu] = min(mrVmin_clu); %center sites

% perform CAR after centering at the center site
miSites_clu = P.miSites(:,viSite_clu);
mrVmin_clu = squeeze_(min(trWav_clu));
[vrVmin_clu, ~] = min(mrVmin_clu);

% cluster specifications
vrVmin_clu = abs(vrVmin_clu);
vrSnr_clu = (vrVmin_clu ./ vrVrms_site(viSite_clu))';
vrThresh_site = vrVrms_site * P.qqFactor;
vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -vrThresh_site(viSite_clu)));

if ~isempty(snr_thresh)
    viClu_keep = find(abs(vrSnr_clu) > snr_thresh);
    [trWav_clu, trWav_raw_clu] = multifun_(@(x)x(:,:,viClu_keep), trWav_clu, trWav_raw_clu);
    [viSite_clu, vrVmin_clu, vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site] = ...
        multifun_(@(x)x(viClu_keep), viSite_clu, vrVmin_clu, vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site);
    vlSpk_keep = ismember(viClu, viClu_keep);
    [S_gt.viClu, S_gt.viTime] = multifun_(@(x)x(vlSpk_keep), S_gt.viClu, S_gt.viTime);
    viMap = 1:max(viClu_keep);
    viMap(viClu_keep) = 1:numel(viClu_keep);
    S_gt.viClu = viMap(S_gt.viClu); % compact, no-gap
end

S_gt = struct_add_(S_gt, trWav_clu, trWav_raw_clu, viSite_clu, vrVmin_clu, ...
    vrSnr_clu, cviSpk_clu, vnSite_clu, vrVrms_site, miSites_clu);
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [via1, via2, via3, via4, vib1, vib2, vib3, vib4] = sgfilt4_(n1, fGpu)
persistent n1_prev_ via1_ via2_ via3_ via4_ vib1_ vib2_ vib3_ vib4_
if nargin<2, fGpu=0; end

% Build filter coeff
if isempty(n1_prev_), n1_prev_ = 0; end
try a = size(via1_); catch, n1_prev_ =0; end
if n1_prev_ ~= n1 %rebuild cache
    vi0 = int32(1:n1);        
    vi0 = gpuArray_(vi0, fGpu);
    via4_ = vi0+4; via4_(end-3:end)=n1;   vib4_ = vi0-4; vib4_(1:4)=1;
    via3_ = vi0+3; via3_(end-2:end)=n1;   vib3_ = vi0-3; vib3_(1:3)=1;
    via2_ = vi0+2; via2_(end-1:end)=n1;   vib2_ = vi0-2; vib2_(1:2)=1;
    via1_ = vi0+1; via1_(end)=n1;         vib1_ = vi0-1; vib1_(1)=1;
    n1_prev_ = n1;
end

% Copy from cache
via4 = via4_;   vib4 = vib4_;
via3 = via3_;   vib3 = vib3_;
via2 = via2_;   vib2 = vib2_;
via1 = via1_;   vib1 = vib1_;
end %func


%--------------------------------------------------------------------------
function [miA, miB, viC] = sgfilt_init_(nData, nFilt, fGpu)
persistent miA_ miB_ viC_ nData_ nFilt_
if nargin<2, fGpu=0; end

% Build filter coeff
if isempty(nData_), nData_ = 0; end
try a = size(miA_); catch, nData_ = 0; end
if nData_ == nData && nFilt_ == nFilt
    [miA, miB, viC] = deal(miA_, miB_, viC_);
else
    vi0 = gpuArray_(int32(1):int32(nData), fGpu)';
    vi1 = int32(1):int32(nFilt);
    miA = min(max(bsxfun(@plus, vi0, vi1),1),nData);
    miB = min(max(bsxfun(@plus, vi0, -vi1),1),nData);
    viC = gpuArray_(int32(-nFilt:nFilt), fGpu);
    [nData_, nFilt_, miA_, miB_, viC_] = deal(nData, nFilt, miA, miB, viC);
end
end %func


%--------------------------------------------------------------------------
function vrVrms_site = mr2rms_(mnWav2, max_sample)
% uses median to estimate RMS
if nargin<2, max_sample = []; end
if isempty(max_sample)
    vrVrms_site = median(abs(mnWav2));
else
    vrVrms_site = median(abs(subsample_mr_(mnWav2, max_sample, 1)));
end
vrVrms_site = single(vrVrms_site) / 0.6745;
end


%--------------------------------------------------------------------------
function [S_clu, S0] = post_merge_(S_clu, P, fPostCluster)
% waveform based merging. find clusters within maxSite
% also removes duplicate spikes
global tnWav_spk tnWav_raw;
if nargin<3, fPostCluster=1; end

nRepeat_merge = get_set_(P, 'nRepeat_merge', 10);
fMerge_pv = 0;
fClean_clu = 1;

if fPostCluster, S_clu = postCluster_(S_clu, P); end

% Add S_clu fields
S0 = get0_();
S_clu = rmfield_(S_clu, 'tmrWav_clu', 'vrPosX_clu', 'vrPosY_clu', 'vrVpp_clu', ...
    'vrVpp_uv_clu', 'vrVmin_uv_clu', 'vrSnr_clu', 'vnSite_clu', ...
    'vrIsoDist_clu', 'vrLRatio_clu', 'vrIsiRatio_clu');
S_clu = S_clu_refresh_(S_clu);
S_clu = S_clu_sort_(S_clu, 'viSite_clu');
S_clu = rmfield_(S_clu, 'csNote_clu');

if isempty(tnWav_raw)
    fprintf(2, 'postmerge: Loading tnWav_raw... (should not see this)\n');
    tnWav_raw = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), 'int16', S0.dimm_raw);
end
if fClean_clu, S_clu = S_clu_cleanup_(S_clu, P); end
S_clu = post_merge_wav_(S_clu, nRepeat_merge, P);

% merge based on pca
if fMerge_pv
    for iRepeat = 1:1
        [S_clu, nMerges_clu] = S_clu_pv_merge_(S_clu, P); %update mrWavCor when you merge
        if nMerges_clu<1, break; end
    end
end
S_clu = S_clu_refresh_(S_clu);
S_clu = S_clu_sort_(S_clu, 'viSite_clu');
S_clu = S_clu_update_wav_(S_clu, S0, P);
                
% set diagonal element
[S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, tnWav_raw));
S_clu.P = P;
S_clu = S_clu_position_(S_clu);
S_clu.csNote_clu = cell(S_clu.nClu, 1);  %reset note
S_clu = S_clu_quality_(S_clu, P);
[S_clu, S0] = S_clu_commit_(S_clu, 'post_merge_');
end %func


%--------------------------------------------------------------------------
function S_clu = post_merge_wav_(S_clu, nRepeat_merge, P)
global tnWav_spk tnWav_raw;

S0 = get0_();
% P = S0.P;

% create covariance matrix (mrDist_wav)
S_clu = clu2wav_(S_clu, S0.viSite_spk, tnWav_spk, tnWav_raw);
S_clu.mrWavCor = S_clu_wavcor_(S_clu, P);  

for iRepeat = 1:nRepeat_merge %single-pass vs dual-pass correction
    [S_clu, nMerges_clu] = S_clu_wavcor_merge_(S_clu, P);
    if nMerges_clu < 1, break; end
end %for
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_update_wav_(S_clu, S0, P)
global tnWav_spk tnWav_raw;
if nargin<3, P = get0_('P'); end
if nargin<2, S0 = get0_(); end
S_clu = clu2wav_(S_clu, S0.viSite_spk, tnWav_spk, tnWav_raw);
S_clu.mrWavCor = S_clu_wavcor_(S_clu, P);  
S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, S_clu_self_corr_(S_clu, tnWav_raw));
end %func


%--------------------------------------------------------------------------
function mr = set_diag_(mr, vr)
n = min(size(mr));
m = numel(vr);
mr(sub2ind([n,n], 1:m, 1:m)) = vr;
end %func


%--------------------------------------------------------------------------
function [vr, vi] = get_diag_(mr)
n = min(size(mr));
vi = sub2ind([n,n], 1:n, 1:n);
vr = mr(vi);
end %func


%--------------------------------------------------------------------------
function batch_(vcFile_batch, vcCommand)
% batch process parameter files (.batch) file
% batch_(myfile.batch, vcCommand): contains a list of .prm files
% batch_(vcFile_batch, vcFile_template): batch contains .bin files

if nargin<2, vcCommand=[]; end
if isempty(vcCommand), vcCommand = 'spikesort'; end

if matchFileExt_(vcCommand, '.prm')
    vcFile_template = vcCommand;
    vcCommand = 'spikesort';
else
    vcFile_template = jrcpath_(read_cfg_('default_prm'));
end

csFiles_prm = load_batch_(vcFile_batch);
vlFile_err = true(size(csFiles_prm));
csErr_file = cell(size(csFiles_prm));

for iFile = 1:numel(csFiles_prm)
    vcFile_prm_ = csFiles_prm{iFile};
    if matchFileExt_(vcFile_prm_, {'.bin', '.dat'}); 
        vcFile_prm_ = makeprm_(vcFile_prm_, vcFile_template, 0);
    end
    try            
        jrc3('clear');
        fprintf('Processing file %d/%d: %s\n', iFile, numel(csFiles_prm), vcFile_prm_); 
        jrc3(vcCommand, vcFile_prm_);
        vlFile_err(iFile) = 0;
    catch
        csErr_file{iFile} = lasterr();
        fprintf(2, 'Error in file %d/%d: %s\n', iFile, numel(csFiles_prm), vcFile_prm_);
        fprintf(2, '\t%s\n\n', csErr_file{iFile});
    end
end %for

% summary
if any(vlFile_err) %error found
    fprintf(2,'Errors found in the .prm files below:\n');
    viFile_err = find(vlFile_err);
    for iErr = 1:numel(viFile_err)
        vcFile_err_ = csFiles_prm{viFile_err(iErr)};
        fprintf(2,'\tError %d/%d: %s\n', iErr, numel(viFile_err), vcFile_err_);
        fprintf(2,'\t\t%s\n\n', csErr_file{iFile});
    end
else
    fprintf('All %d files processed successfully\n', numel(csFiles_prm));
end
edit(vcFile_batch);
end %func


%--------------------------------------------------------------------------
function [mnWav1, vnWav1_mean] = wav_car_(mnWav1, P)
% take common average referencing (CAR) on the filtered trace (mnWav1)
% fprintf('Common average referencing (CAR)\n\t'); t1=tic;
fRepairSites = 0; % bad sites get repaired by averaging vertical neighbors
vnWav1_mean = [];
if strcmpi(P.vcCommonRef, 'tmean') || strcmpi(P.vcCommonRef, 'nmean')
    trimLim = [.25, .75];    
    maxSite_ref = (P.nSites_ref + P.nSites_excl_ref - 1)/2;
    miSite_ref = findNearSites_(P.mrSiteXY, maxSite_ref, P.viSiteZero, P.viShank_site);
    miSite_ref = miSite_ref(P.nSites_excl_ref+1:end, :); %excl three nearest sites
    viChan_keep = round(trimLim * size(miSite_ref,1));
    viChan_keep = (viChan_keep(1)+1):viChan_keep(2);
    mnWav1_pre = mnWav1;
    if strcmpi(P.vcCommonRef, 'tmean')
        for iChan=1:size(mnWav1,2)
            mnWav2 = sort(mnWav1_pre(:, miSite_ref(:,iChan)), 2);
            gvr_tmean = sum(mnWav2(:,viChan_keep), 2); %may go out of range
            gvr_tmean = int16(single(gvr_tmean)/numel(viChan_keep));
            mnWav1(:,iChan) = mnWav1_pre(:,iChan) - gvr_tmean;
            fprintf('.');
        end
    else
        for iChan=1:size(mnWav1,2)
            gvr_tmean = sum(mnWav1_pre(:, miSite_ref(:,iChan)), 2); %may go out of range
            gvr_tmean = int16(single(gvr_tmean)/size(miSite_ref,1));
            mnWav1(:,iChan) = mnWav1_pre(:,iChan) - gvr_tmean;
            fprintf('.');
        end
    end
elseif strcmpi(P.vcCommonRef, 'mean') %compute mean excl. viSiteZero    
    vnWav1_mean = mean_excl_(mnWav1, P);
    mnWav1 = bsxfun(@minus, mnWav1, vnWav1_mean);
end
% viSiteZero should be treated carefully. try to repair using nearest sites?
if get_(P, 'fMeanSite_drift')
    mnWav1 = meanSite_drift_(mnWav1, P); 
elseif fRepairSites
    mnWav1 = meanSite_drift_(mnWav1, P, P.viSiteZero); 
else
    mnWav1(:,P.viSiteZero) = 0;
end
% fprintf('\n\ttook %0.1fs.\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function batch_verify_(vcFile_batch, vcCommand)
% batch process parameter files (.batch) file
% Example
%   jrc3 batch-verify skip my.batch
%       just does the verification plot for all files in .batch file
if ~exist(vcFile_batch, 'file'), fprintf(2, 'File does not exist\n'); return; end
edit(vcFile_batch); %show script
if nargin<2, vcCommand=[]; end
if isempty(vcCommand), vcCommand='spikesort'; end
csFiles_prm = load_batch_(vcFile_batch);

if ~strcmpi(vcCommand, 'skip')
    for iFile=1:numel(csFiles_prm)
        try
            vcFile_prm1 = csFiles_prm{iFile};
            jrc3('clear');
            jrc3(vcCommand, vcFile_prm1);
            if isempty(strfind(vcCommand, 'verify'))
                jrc3('verify', vcFile_prm1); % try silent verify and collect result
            end
        catch
            disperr_();
        end
    end %for
end
fprintf('\nSummary for %s\n', vcFile_batch);
% Collect data
[csSnr, csFp, csFn, cvnSite] = deal(cell(size(csFiles_prm)));
for iFile=1:numel(csFiles_prm)    
    try
        S_score1 = load(strrep(csFiles_prm{iFile}, '.prm', '_score.mat'));    
        csSnr{iFile} = gather_(S_score1.vrSnr_min_gt');    
        csFp{iFile} = S_score1.S_score_clu.vrFp;
        csFn{iFile} = S_score1.S_score_clu.vrMiss;        
        disp(csFiles_prm{iFile});
        disp_score_(csSnr{iFile}, csFp{iFile}, csFn{iFile});
        cvnSite{iFile} = S_score1.vnSite_gt;
    catch
        disperr_();
    end
end

[vrSnr, vrFp, vrFn, vnSite] = multifun_(@(x)cell2mat_(x'), csSnr, csFp, csFn, cvnSite);
disp('All files pooled:');
disp_score_(vrSnr, vrFp, vrFn);

% Plot
vpp_bin = 2; %25 for vpp

figure;  ax=[];
ax(1)=subplot(311); 
boxplot_(vrFp(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Positive'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
ylim([0, .2]); xlim([0 40]);
title_(vcFile_batch);

ax(2)=subplot(312); 
boxplot_(vrFn(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('False Negative'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
ylim([0, .2]); xlim([0 40]);
set(gcf,'Color','w');

ax(3)=subplot(313); 
boxplot_(vnSite(:), vrSnr(:), vpp_bin, [3, 40]); ylabel('#sites>thresh'); 
xlabel('SNR (Vp/Vrms)'); grid on; set(gca,'YScale','linear');
linkaxes(ax,'x');
set(gcf,'Color','w');
ylim([0, 16]); 
end %func


%--------------------------------------------------------------------------
function disp_score_(vrSnr, vrFp, vrFn, snr_thresh)
if nargin<4, snr_thresh = 10; end %use 10 as a default
% disp_score_(vrSnr, vrFp, vrFn, snr_thresh)
fprintf('SNR>%d Groundtruth Units\n', snr_thresh);
viSnr = find(vrSnr > snr_thresh);
fprintf('\tfalse positive (%%): '); disp_stats_(vrFp(viSnr)*100);
fprintf('\tfalse negative (%%): '); disp_stats_(vrFn(viSnr)*100);
fprintf('\tscore (1-FP-FN) (%%): '); disp_stats_((1-vrFp(viSnr)-vrFn(viSnr))*100);
end %func


%--------------------------------------------------------------------------
function import_gt_silico_(vcFile_mat)
if matchFileExt_(vcFile_mat, '.mat')
    % catalin's raster function
    S = load(vcFile_mat);    
    vnSpk = cellfun(@numel, S.a);    
    viClu = int32(cell2mat_(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0)));
    viTime = int32(cell2mat_(S.a) * 20); % Convert to sample # (saved in ms unit & sampling rate =20KHZ)
end
S_gt = makeStruct_(viClu, viTime);
vcFile_gt = subsFileExt_(vcFile_mat, '_gt.mat');
write_struct_(vcFile_gt, S_gt);
end %func


%--------------------------------------------------------------------------
function close_hFig_traces_(hFig, event)
try
    if ~ishandle(hFig), return; end
    if ~isvalid(hFig), return; end
    S_fig = get(hFig, 'UserData');    
    fclose_(S_fig.fid_bin);
    try delete(hFig); catch; end %close one more time
catch
    disperr_();
    close(hFig);
end
end %func


%--------------------------------------------------------------------------
function fclose_(fid, fVerbose)
if nargin<2, fVerbose = 0; end
if isempty(fid), return; end
if ischar(fid), return; end
try 
    fclose(fid); 
    if fVerbose, disp('File closed.'); end
catch
    disperr_(); 
end
end %func


%--------------------------------------------------------------------------
function hTitle = title_(hAx, vc)
% title_(vc)
% title_(hAx, vc)

if nargin==1, vc=hAx; hAx=[]; end
% Set figure title

if isempty(hAx), hAx = gca; end
hTitle = get_(hAx, 'Title');
if isempty(hTitle)
    hTitle = title(hAx, vc, 'Interpreter', 'none', 'FontWeight', 'normal');
else
    set_(hTitle, 'String', vc, 'Interpreter', 'none', 'FontWeight', 'normal');
end
end %func


%--------------------------------------------------------------------------
function [mnWav2, vnWav2_mean] = filt_car_(mnWav2, P, mnWav1_pre, mnWav1_post, fTrim_pad)
% Apply filter and CAR
% @TODO: edge case
if nargin<3, mnWav1_pre = []; end
if nargin<4, mnWav1_post = []; end
if nargin<5, fTrim_pad = 1; end
n_pre = size(mnWav1_pre,1);
n_post = size(mnWav1_post,1);
if n_pre > 0 || n_post > 0
    mnWav2 = [mnWav1_pre; mnWav2; mnWav1_post];
end
P.vcFilter = get_filter_(P);
switch lower(P.vcFilter)
    case {'sgdiff', 'sgfilt'}
        mnWav2 = sgfilt_(mnWav2, P.nDiff_filt, P.fGpu);
    case 'bandpass'
        try
            mnWav2 = filtfilt_chain(single(mnWav2), P);
        catch
            fprintf('GPU filtering failed. Trying CPU filtering.\n');
            mnWav2 = filtfilt_chain(single(mnWav2), setfield(P, 'fGpu', 0));
        end
        mnWav2 = int16(mnWav2);
    otherwise
        error('filt_car_: invalid filter option (vcFilter=''%s'')', P.vcFilter);
end  %switch

% trim padding
if (n_pre > 0 || n_post > 0) && fTrim_pad
    mnWav2 = mnWav2(n_pre+1:end-n_post,:);
end
% mean subtract
% if strcmpi(P.vcCommonRef, 'mean')
%     mnWav2 = bsxfun(@minus, mnWav2, int16(mean(mnWav2,2)));
% end
% % Better to do CAR after spike detection and merging
% if ~strcmpi(P.vcCommonRef, 'mean')
[mnWav2, vnWav2_mean] = wav_car_(mnWav2, P); %global subtraction before 
% end
end %func


%--------------------------------------------------------------------------
function mnWav1 = fread_(fid_bin, dimm_wav, vcDataType)
% Get around fread bug (matlab) where built-in fread resize doesn't work
try
    if isempty(dimm_wav)
        mnWav1 = fread(fid_bin, inf, ['*', vcDataType]);
    else
        mnWav1 = fread(fid_bin, prod(dimm_wav), ['*', vcDataType]);
        if numel(mnWav1) == prod(dimm_wav)
            mnWav1 = reshape(mnWav1, dimm_wav);
        else
            dimm2 = floor(numel(mnWav1) / dimm_wav(1));
            if dimm2 >= 1
                mnWav1 = reshape(mnWav1, dimm_wav(1), dimm2);
            else
                mnWav1 = [];
            end
        end
    end
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function frewind_(fid_bin, dimm_wav, vcDataType)
% move the file pointer back by the dimm_wav, vcDatatype
fseek(fid_bin, -1 * prod(dimm_wav) * bytesPerSample_(vcDataType), 'cof');
end %func


%--------------------------------------------------------------------------
function csDesc = describe_(vcFile_prm)
%describe_()
%describe_(vcFile_prm)
%describe_(S0)

% describe _jrc.mat file
if nargin==0
    S0 = get(0, 'UserData'); 
elseif isstruct(vcFile_prm)
    S0 = vcFile_prm;
else
    S0 = load(strrep(vcFile_prm, '.prm', '_jrc.mat'));
end
P = S0.P;

nSites = numel(P.viSite2Chan);
tDur = double(max(S0.viTime_spk) - min(S0.viTime_spk)) / P.sRateHz;
nSpk = numel(S0.viTime_spk);
nSitesPerEvent = P.maxSite*2+1;


csDesc = {};
    csDesc{end+1} = sprintf('Recording file');
    csDesc{end+1} = sprintf('    Recording file          %s', P.vcFile);
    csDesc{end+1} = sprintf('    Probe file              %s', P.probe_file);
    csDesc{end+1} = sprintf('    Recording Duration      %0.1fs ', tDur);
    csDesc{end+1} = sprintf('    #Sites                  %d', nSites);
    csDesc{end+1} = sprintf('Events');
    csDesc{end+1} = sprintf('    #Spikes                 %d', nSpk);
    csDesc{end+1} = sprintf('    Feature                 %s', P.vcFet);
    csDesc{end+1} = sprintf('    #Sites/event            %d', nSitesPerEvent);
if ~isempty(get_(S0, 'nLoads'))
    csDesc{end+1} = sprintf('    #Loads                  %d', S0.nLoads);
end

if isfield(S0, 'S_clu')
    S_clu = S0.S_clu;
    csDesc{end+1} = sprintf('Cluster');    
    csDesc{end+1} = sprintf('    #Clusters               %d', S_clu.nClu);
    csDesc{end+1} = sprintf('    #Unique events          %d', sum(S_clu.viClu>0));
    csDesc{end+1} = sprintf('    min. spk/clu            %d', P.min_count);
    csDesc{end+1} = sprintf('    Cluster run-time        %0.1fs', S_clu.t_runtime);    
end
try
    runtime_total = S0.runtime_detect + S0.runtime_sort;
    csDesc{end+1} = sprintf('Runtime (s)');
    csDesc{end+1} = sprintf('    Detect + feature        %0.1fs', S0.runtime_detect);    
    csDesc{end+1} = sprintf('    Cluster                 %0.1fs', S0.runtime_sort);    
    csDesc{end+1} = sprintf('    Total                   %0.1fs', runtime_total);
    csDesc{end+1} = sprintf('    Runtime speed           x%0.1f realtime', tDur / runtime_total);
catch
    ;
end

if nargout==0
    cellfun(@(x)disp(x), csDesc);
end

end %func


%--------------------------------------------------------------------------
function [viTime1, viSite1, viSpk1] = trimSpikes_(nlim_show)
error('not implemented');

viTime1=[];  viSite1=[];
if ~isfield(P, 'viSpk') || ~isfield(P, 'viSite'), return; end
ilim = round(P.tlim * P.sRateHz);
viSpk1 = find(P.viSpk >= ilim(1) & P.viSpk < ilim(end));
viTime1 = P.viSpk(viSpk1) - ilim(1);
viSite1 = P.viSite(viSpk1);
end %func


%--------------------------------------------------------------------------
function delete_multi_(varargin)
% provide cell or multiple arguments
for i=1:nargin
    try
        vr1 = varargin{i};
        if numel(vr1)==1
            delete(varargin{i}); 
        elseif iscell(vr1)
            for i1=1:numel(vr1)
                try
                    delete(vr1{i1});
                catch
                end
            end
        else
            for i1=1:numel(vr1)
                try
                    delete(vr1(i1));
                catch
                end
            end
        end
    catch
    end
end
end %func


%--------------------------------------------------------------------------
function vi = cell2vi_(cvi)
% convert cell index to array of index
vn_site = cellfun(@(x)numel(x), cvi); %create uniform output
vi = cell(numel(cvi), 1);
for iSite=1:numel(cvi)
    vi{iSite} = iSite * ones(vn_site(iSite), 1);
end
vi = cell2mat_(vi);
end %func


%--------------------------------------------------------------------------
function out = inputdlg_num_(vcGuide, vcVal, vrVal)
csAns = inputdlg_(vcGuide, vcVal, 1, {num2str(vrVal)});
try
    out = str2double(csAns{1});
catch
    out = nan;
end
end %func


%--------------------------------------------------------------------------
function manual_(P, vcMode, vcFile_ui)
% display manual sorting interface
% global tnWav_spk tnWav_raw mrFet; 
global fDebug_ui

if nargin<2, vcMode = 'normal'; end %{'normal', 'debug'}
if nargin<3, vcFile_ui = ''; end

% Load info
S0 = load_cached_(P);
fDebug_ui = 0;
P.fGpu = 0; %do not use GPU for manual use
set0_(fDebug_ui, P);
switch lower(vcMode)
    case 'normal'
        switch lower(questdlg_('Load last saved?', 'Confirmation'))
            case 'no'
                [S_clu, S0] = post_merge_(S0.S_clu, P);
                S0 = clear_log_(S0);
            case 'cancel'
                return;
            case 'yes'
                S0 = set0_(P); %update the P structure
                S0.S_clu = S_clu_update_wav_(S0.S_clu, S0, P);                
        end
        
    case 'debug'
        fDebug_ui = 1;
        S0 = set0_(fDebug_ui);
        [S_clu, S0] = post_merge_(S0.S_clu, P); %redo the clustering (reset to auto)
        S0 = set0_(P);
        
    case 'record'
        % Generate a script file and save it to mat file
        % fully self-contained. debug purpose only. Contains _jrc + exta
        % name it as _gui.mat
        fprintf(2, 'manual-record not implemented.\n'); return;
    case 'replay'
        % replay the recorded demo. Need a script file
        fprintf(2, 'manual-reply not implemented.\n'); return;
end

% Create figures
hMsg = msgbox_('Plotting...'); t1=tic;
set(0, 'UserData', S0);
S0 = figures_manual_(P); %create figures for manual interface
clear mouse_figure;
clear get_fig_cache_; %clear persistent figure handles

% Set fields
S0 = struct_merge_(S0, ...
    struct('iCluCopy', 1, 'iCluPaste', [], 'hCopy', [], 'hPaste', [], 'nSites', numel(P.viSite2Chan)));
set(0, 'UserData', S0);

% hFigRD
S0.S_clu = plot_FigRD_(S0.S_clu, P); % ask user before doing so

% Set initial amplitudes
set(0, 'UserData', S0);
plot_FigWavCor_(S0);  % hFigWavCor
S0 = plot_FigWav_(S0); % hFigWav %do this after for ordering

% hFigProj, hFigHist, hFigIsi, hFigCorr, hFigPos, hFigMap, hFigTime
close_(get_fig_('FigTrial')); %close previous FigTrial figure
S0 = button_CluWav_simulate_(1, [], S0); %select first clu
auto_scale_proj_time_();
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom
%S0.cS_log = load_(strrep(P.vcFile_prm, '.prm', '_log.mat'), 'cS_log', 0); 
S_log = load_(strrep(P.vcFile_prm, '.prm', '_log.mat'), [], 0);
if ~isempty(S_log), S0.cS_log = {S_log}; end
S0 = save_log_('start', S0); %crash proof log
set(0, 'UserData', S0);

% Finish up
close_(hMsg);
fprintf('UI creation took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function S0 = figures_manual_(P)
% 'iFig', [], 'name', '', 'pos', [], 'fToolbar', 0, 'fMenubar', 0);
create_figure_('FigPos', [0 0 .15 .5], ['Unit position; ', P.vcFile], 1, 1);
create_figure_('FigMap', [0 .5 .15 .5], ['Probe map; ', P.vcFile], 1, 1);  

create_figure_('FigWav', [.15 .2 .35 .8],['Averaged waveform: ', P.vcFile], 0, 1);
create_figure_('FigTime', [.15 0 .7 .2], ['Time vs. Amplitude; (Sft)[Up/Down] channel; [h]elp; [a]uto scale; ', P.vcFile]);

create_figure_('FigProj', [.5 .2 .35 .5], ['Amplitude projection: ', P.vcFile]);
create_figure_('FigWavCor', [.5 .7 .35 .3], ['Waveform correlation (click): ', P.vcFile]);

create_figure_('FigHist', [.85 .75 .15 .25], ['ISI Histogram: ', P.vcFile]);
create_figure_('FigIsi', [.85 .5 .15 .25], ['Return map: ', P.vcFile]);
create_figure_('FigCorr', [.85 .25 .15 .25], ['Time correlation: ', P.vcFile]);
create_figure_('FigRD', [.85 0 .15 .25], ['Cluster rho-delta: ', P.vcFile]);

csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
cvrFigPos0 = cellfun(@(vc)get(get_fig_(vc), 'OuterPosition'), csFig, 'UniformOutput', 0);
S0 = set0_(cvrFigPos0, csFig);
end %func


%--------------------------------------------------------------------------
function S0 = button_CluWav_simulate_(iCluCopy, iCluPaste, S0)
if nargin<3,  S0 = get(0, 'UserData'); end
if nargin<2, iCluPaste = []; end
if iCluCopy == iCluPaste, iCluPaste = []; end
hFig = gcf;
figure_wait_(1, hFig);

S0 = update_cursor_(S0, iCluCopy, 0);
S0 = update_cursor_(S0, iCluPaste, 1);
    
% set(0, 'UserData', S0);
S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'c','t','j','i','v','e','f'}, S0); %'z' to recenter
set(0, 'UserData', S0);

% update PSTH if displayed
plot_raster_(); %psth
auto_scale_proj_time_(); %autoscale
% drawnow;
figure_wait_(0, hFig);
end


%--------------------------------------------------------------------------
function S0 = update_cursor_(S0, iClu, fPaste)
if isempty(iClu), return; end
if isempty(S0), S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
% [hFig, S_fig] = get_fig_cache_('FigWav');

if ~isfield(S0, 'hCopy'), S0.hCopy = []; end
if ~isfield(S0, 'hPaste'), S0.hPaste = []; end

if ~fPaste
    iCluCopy = iClu;    
    if iCluCopy <1 || iCluCopy > S_clu.nClu, return; end
    update_plot_(S0.hPaste, nan, nan); %hide paste
    S0.iCluPaste = []; 
    [S0.iCluCopy, S0.hCopy] = plot_tmrWav_clu_(S0, iCluCopy, S0.hCopy, [0 0 0]);
else
    iCluPaste = iClu;    
    if iCluPaste < 1 || iCluPaste > S_clu.nClu || S0.iCluCopy == iCluPaste, return; end
    [S0.iCluPaste, S0.hPaste] = plot_tmrWav_clu_(S0, iCluPaste, S0.hPaste, [1 0 0]);
end
% set(hFig, 'UserData', S_fig);
cursor_FigWavCor_(S0);
if nargout==0, set(0, 'UserData', S0); end
end %func


%--------------------------------------------------------------------------
function cursor_FigWavCor_(S0)
if nargin==0, S0 = get(0, 'UseData'); end
P = S0.P; S_clu = S0.S_clu;

[hFig, S_fig] = get_fig_cache_('FigWavCor');
if isempty(S_fig)
    [hFig, S_fig] = plot_FigWavCor_(S0);
end
iClu1 = S0.iCluCopy;
if isempty(S0.iCluPaste)
    iClu2 = S0.iCluCopy;
else
    iClu2 = S0.iCluPaste;
end

cor12 = S_clu.mrWavCor(iClu1, iClu2);
set(S_fig.hCursorV, 'XData', iClu1*[1,1], 'YData', [.5, S_clu.nClu+.5]);
title_(S_fig.hAx, sprintf('Clu%d vs. Clu%d: %0.3f; %s', iClu1, iClu2, cor12, S_fig.vcTitle));    
if iClu1==iClu2, color_H = [0 0 0]; else color_H = [1 0 0]; end
set(S_fig.hCursorH, 'YData', iClu2*[1,1], 'XData', [.5, S_clu.nClu+.5], 'Color', color_H);
xlim(S_fig.hAx, trim_lim_(iClu1 + [-6,6], [.5, S_clu.nClu+.5]));
ylim(S_fig.hAx, trim_lim_(iClu2 + [-6,6], [.5, S_clu.nClu+.5]));
end %func


%--------------------------------------------------------------------------
function nFailed = unit_test_(vcArg1, vcArg2, vcArg3)
% 2017/2/24. James Jun. built-in unit test suite (request from Karel Svoboda)
% run unit test
%[Usage]
% unit_test()
%   run all
% unit_test(iTest)
%   run specific test again and show profile
% unit_test('show')
%   run specific test again and show profile
% @TODO: test using multiple datasets and parameters.
global fDebug_ui;

if nargin<1, vcArg1 = ''; end
if nargin<2, vcArg2 = ''; end
if nargin<3, vcArg3 = ''; end

if ~exist_file_('sample.bin'), jrc3('download', 'sample'); end

nFailed = 0;
profile('clear'); %reset profile stats
csCmd = {...  
    'close all; clear all;', ... %start from blank
    'jrc3 clear sample_sample.prm', ...    
    'jrc3 compile', ...
    'jrc3 probe sample.prb', ...
    'jrc3 makeprm sample_list.txt sample.prb', ...
    'jrc3 makeprm sample.bin sample.prb', ...
    'jrc3 probe sample_sample.prm', ...    
    'jrc3 traces sample_sample.prm', 'jrc3 traces', ...
    'jrc3 traces-test sample_sample.prm', ...    
    'jrc3 preview-test sample_sample.prm', ...    
    'jrc3 makeprm sample.bin sample.prb', ...
    'jrc3 detectsort sample_sample.prm', 'jrc3 clear sample_sample.prm', ...
    'jrc3 probe sample_sample_merge.prm', 'jrc3 detectsort sample_sample_merge.prm', 'jrc3 clear sample_sample_merge.prm', ... %multishank, multifile test
    'jrc3 detect sample_sample.prm', 'jrc3 sort sample_sample.prm', ...        
    'jrc3 export-csv sample_sample.prm', ...
    'jrc3 export-spkwav sample_sample.prm', ...
    'jrc3 export-spkwav sample_sample.prm 1', ...
    'jrc3 export-spkamp sample_sample.prm', ...
    'jrc3 export-spkamp sample_sample.prm 1', ...
    'jrc3 export-jrc1 sample_sample.prm', ...
    'jrc3 export-fet sample_sample.prm', ...
    'jrc3 plot-activity sample_sample.prm', ... %     'jrc3 kilosort sample_sample.prm', ...
    'jrc3 clear sample_sample.prm', ...
    'jrc3 traces-test sample_sample.prm', ...
    'jrc3 detectsort sample_sample.prm', ...    
    'jrc3 auto sample_sample.prm', ...    
    'jrc3 manual-test sample_sample.prm', ...
    }; %last one should be the manual test

if ~isempty(vcArg1)
    switch lower(vcArg1)
        case {'show', 'info', 'list', 'help'}
            arrayfun(@(i)fprintf('%d: %s\n', i, csCmd{i}), 1:numel(csCmd)); 
            return;
        case {'manual', 'ui', 'ui-manual'}
            iTest = numel(csCmd); % + [-1,0];
        case {'traces', 'ui-traces'}
            iTest = numel(csCmd)-2; % second last
        otherwise
            iTest = str2num(vcArg1);
    end        
    fprintf('Running test %s: %s\n', vcArg1, csCmd{iTest});
    csCmd = csCmd(iTest);
end

vlPass = false(size(csCmd));
[csError, cS_prof] = deal(cell(size(csCmd)));
vrRunTime = zeros(size(csCmd));
for iCmd = 1:numel(csCmd)   
    eval('close all; fprintf(''\n\n'');'); %clear memory
    fprintf('Test %d/%d: %s\n', iCmd, numel(csCmd), csCmd{iCmd}); 
    t1 = tic;        
    profile('on');
    fDebug_ui = 1;
    set0_(fDebug_ui);
    try
        if any(csCmd{iCmd} == '(' | csCmd{iCmd} == ';') %it's a function        
            evalin('base', csCmd{iCmd}); %run profiler 
        else % captured by profile
            csCmd1 = strsplit(csCmd{iCmd}, ' ');
            feval(csCmd1{:});
        end
        vlPass(iCmd) = 1; %passed test        
    catch
        csError{iCmd} = lasterr();
        fprintf(2, '\tTest %d/%d failed\n', iCmd, numel(csCmd));
    end
    vrRunTime(iCmd) = toc(t1);
    cS_prof{iCmd} = profile('info');    
end
nFailed = sum(~vlPass);

fprintf('Unit test summary: %d/%d failed.\n', sum(~vlPass), numel(vlPass));
for iCmd = 1:numel(csCmd)
    if vlPass(iCmd)
        fprintf('\tTest %d/%d (''%s'') took %0.1fs.\n', iCmd, numel(csCmd), csCmd{iCmd}, vrRunTime(iCmd));
    else
        fprintf(2, '\tTest %d/%d (''%s'') failed:%s\n', iCmd, numel(csCmd), csCmd{iCmd}, csError{iCmd});
    end        
end

if numel(cS_prof)>1
    assignWorkspace_(cS_prof);
    disp('To view profile, run: profview(0, cS_prof{iTest});');
else
    profview(0, cS_prof{1});
end
fDebug_ui = [];
set0_(fDebug_ui);
end %func


%--------------------------------------------------------------------------
function S_clu = plot_FigRD_(S_clu, P)
% P = funcDefStr_(P, ...
%     'delta1_cut', .5, 'rho_cut', -2.5, 'fDetrend', 0, 'fAskUser', 1, ...
%     'min_count', 50, 'y_max', 2, 'rhoNoise', 0, 'minGamma', -1, 'fNormDelta', 0, 'fExclMaxRho', 0, 'fLabelClu', 1);
% P.y_max = 1;
% P.fAskUser = 1;

[hFig, S_fig] = get_fig_cache_('FigRD');
figure(hFig); clf;

if isfield(S_clu, 'cS_clu_shank')
    cellfun(@(S_clu1)plot_FigRD_(S_clu1, P), S_clu.cS_clu_shank);
    return;
end

if isempty(P.delta1_cut), P.delta1_cut = S_clu.P.delta1_cut; end
if isempty(P.rho_cut), P.rho_cut = S_clu.P.rho_cut; end
if isempty(P.min_count), P.min_count = S_clu.P.min_count; end
if ~isfield(P, 'vcDetrend_postclu'), P.vcDetrend_postclu = 'none'; end

switch P.vcDetrend_postclu
    case 'none'
        icl = find(S_clu.rho(:) > 10^(P.rho_cut) & S_clu.delta(:) > 10^(P.delta1_cut));
        x = log10(S_clu.rho(:));
        y = log10(S_clu.delta(:));
        fDetrend = 0;
    case 'global'
%         [icl, x, y] = detrend_ztran_(S_clu.rho, S_clu.delta, P.rho_cut, P.delta1_cut);
        [icl, x, y] = detrend_local_(S_clu, P, 0);
%         y(y<=0) = nan;
%         y = log10(y);
%         y(~isreal(y)) = nan;
        y = log10(y);
        y(~isreal(y)) = nan;
        fDetrend = 1;
    case 'local'
%         icl = S_clu.icl;
%         x = log10(S_clu.rho(:));
%         y = S_clu.delta(:);
%         viDetrend = find(S_clu.delta < 1 & S_clu.rho > 10^P.rho_cut & S_clu.rho < .1 & isfinite(x) & isfinite(y));
%         [~,y] = detrend_(x,y, viDetrend);
        [icl, x, y] = detrend_local_(S_clu, P, 1);
        y = log10(y);
        y(~isreal(y)) = nan;
        fDetrend = 1;
end

hold on; plot(x, y, '.');
axis tight;
ylim([-.5, 2]);
set(gcf,'color','w');
set(gcf, 'UserData', struct('x', x, 'y', y)); grid on; 
set(gca,'XScale','linear', 'YScale', 'linear');
plot(P.rho_cut*[1 1], get(gca,'YLim'), 'r--', get(gca,'XLim'), P.delta1_cut*[1, 1], 'r--');
xlabel('log10 rho'); ylabel(sprintf('log10 delta (detrend=%d)', fDetrend));

% label clusters
if isfield(S_clu, 'icl')
    icl = S_clu.icl; % do not overwrite
end
x_icl = double(x(icl));
y_icl = double(y(icl));
% if P.fLabelClu
%     arrayfun(@(i)text(x_icl(i), y_icl(i), sprintf('%dn%d',i,S_clu.vnSpk_clu(i)), 'VerticalAlignment', 'bottom'), 1:numel(icl));
% end
hold on; 
plot(x_icl, y_icl, 'r.');
grid on; 
nClu = numel(unique(S_clu.viClu(S_clu.viClu>0))); %numel(icl)
title_(sprintf('#clu:%d, rho-cut:%f, delta-cut:%f', nClu, P.rho_cut, P.delta1_cut));
drawnow;

end %func


%--------------------------------------------------------------------------
function S0 = plot_FigWav_(S0)
if nargin<1, S0 = get0_(); end
P = S0.P; S_clu = S0.S_clu;

[hFig, S_fig] = get_fig_cache_('FigWav'); 

% Show number of spikes per clusters
% hold on; tight_plot(gca, [.04 .04], [.04 .02]);
P.LineWidth = 1; %plot a thicker line
P.viSite_clu = S_clu.viSite_clu;
nSites = numel(P.viSite2Chan);
if isempty(S_fig)
    % initialize
    S_fig.maxAmp = P.maxAmp;
    S_fig.hAx = axes_new_(hFig);
    set(gca, 'Position', [.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual'); 
    xlabel('Cluster #');    ylabel('Site #');   grid on;    
    S_fig.vcTitle = 'Scale: %0.1f uV; [H]elp; [Left/Right]:Select cluster; (Sft)[Up/Down]:scale; [M]erge; [S]plit auto; [D]elete; [A]:Resample spikes; [P]STH; [Z]oom; in[F]o; [Space]:Find similar';
    title_(sprintf(S_fig.vcTitle, S_fig.maxAmp)); %update scale
    
%     set(gca, 'ButtonDownFcn', @(src,event)button_CluWav_(src,event), 'BusyAction', 'cancel');
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWav_, 'CloseRequestFcn', @exit_manual_, 'BusyAction', 'cancel');
    axis([0, S_clu.nClu + 1, 0, nSites + 1]);
    add_menu_(hFig, P);      
    mouse_figure(hFig, S_fig.hAx, @button_CluWav_);
    S_fig = plot_spkwav_(S_fig, P); %plot spikes
    S_fig = plot_tnWav_clu_(S_fig, P); %do this after plotSpk_
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hSpkAll, S_fig);
else
%     mh_info = [];
    S_fig = plot_spkwav_(S_fig, P); %plot spikes
    try delete(S_fig.vhPlot); catch; end %delete old text
    S_fig = rmfield_(S_fig, 'vhPlot');
    S_fig = plot_tnWav_clu_(S_fig, P); %do this after plotSpk_
end

% create text
% S0 = set0_(mh_info);
S_fig = figWav_clu_count_(S_fig);
S_fig.csHelp = { ...            
    '[Left-click] Cluter select/unselect (point at blank)', ...
    '[Right-click] Second cluster select (point at blank)', ...
    '[Pan] hold wheel and drag', ...
    '[Zoom] mouse wheel', ...
    '[X + wheel] x-zoom select', ...
    '[Y + wheel] y-zoom select', ...
    '[SPACE] clear zoom', ...
    '[(shift) UP]: increase amplitude scale', ...
    '[(shift) DOWN]: decrease amplitude scale', ...
    '------------------', ...
    '[H] Help', ...       
    '[S] Split auto', ...
    '[W] Spike waveforms (toggle)', ...                                  
    '[M] merge cluster', ...
    '[D] delete cluster', ...
    '[A] Resample spikes', ...
    '[N] Spike count (toggle)', ...
    '[Z] zoom selected cluster', ...
    '[R] reset view', ...
    '------------------', ...
    '[U] update all', ...  
    '[C] correlation plot', ...            
    '[T] show amp drift vs time', ...            
    '[J] projection view', ...            
    '[V] ISI return map', ...
    '[I] ISI histogram', ...
    '[E] Intensity map', ...
    '[P] PSTH display', ...
    '[O] Overlap average waveforms across sites', ...
    }; 
set(hFig, 'UserData', S_fig);
xlabel('Clu #'); ylabel('Site #');
end %func


%--------------------------------------------------------------------------
function add_menu_(hFig, P)
drawnow;
posvec = get(hFig, 'OuterPosition');

set(hFig, 'MenuBar','None');
mh_file = uimenu(hFig,'Label','File'); 
uimenu(mh_file,'Label', 'Save', 'Callback', @save_manual_);
uimenu(mh_file,'Label', 'Save figures as .fig', 'Callback', @(h,e)save_figures_('.fig'));
uimenu(mh_file,'Label', 'Save figures as .png', 'Callback', @(h,e)save_figures_('.png'));
uimenu(mh_file,'Label', 'Describe', 'Callback', @(h,e)msgbox_(describe_()));
uimenu(mh_file,'Label', 'Edit prm file', 'Callback', @edit_prm_);
uimenu(mh_file,'Label', 'Reload prm file', 'Callback', @reload_prm_);
uimenu(mh_file,'Label', 'Export units to csv', 'Callback', @export_csv_);
uimenu(mh_file,'Label', 'Export all mean unit waveforms', 'Callback', @export_tmrWav_clu_);
uimenu(mh_file,'Label', 'Export selected mean unit waveforms', 'Callback', @(h,e)export_mrWav_clu_);
uimenu(mh_file,'Label', 'Export all waveforms from the selected unit', 'Callback', @(h,e)export_tnWav_spk_);
uimenu(mh_file,'Label', 'Export firing rate for all units', 'Callback', @(h,e)export_rate_);
uimenu(mh_file,'Label', 'Exit', 'Callback', @exit_manual_);

mh_edit = uimenu(hFig,'Label','Edit'); 
uimenu(mh_edit,'Label', '[M]erge', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'm'));
uimenu(mh_edit,'Label', '[D]elete', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'd'));
uimenu(mh_edit,'Label', 'Delete auto', 'Callback', @(h,e)delete_auto_());
uimenu(mh_edit,'Label', '[S]plit', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'));
uimenu(mh_edit,'Label', 'Auto split max-chan', 'Callback', @(h,e)auto_split_(0));
uimenu(mh_edit,'Label', 'Auto split multi-chan', 'Callback', @(h,e)auto_split_(1));
% uimenu(mh_edit,'Label', 'Restore last deleted', 'Callback', @(h,e)restore_clu_());
uimenu(mh_edit,'Label', 'Annotate', 'Callback', @(h,e)unit_annotate_());

mh_view = uimenu(hFig,'Label','View'); 
uimenu(mh_view,'Label', 'Show traces', 'Callback', @(h,e)traces_());
uimenu(mh_view,'Label', 'View all [R]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'r'));
uimenu(mh_view,'Label', '[Z]oom selected', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'z'));
uimenu(mh_view,'Label', '[W]aveform (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'w'));
uimenu(mh_view,'Label', '[N]umbers (toggle)', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'n'));
uimenu(mh_view,'Label', 'Show raw waveform', 'Callback', @(h,e)raw_waveform_(h), ...
    'Checked', ifeq_(get_(P, 'fWav_raw_show'), 'on', 'off'));
%uimenu(mh_view,'Label', 'Threshold by sites', 'Callback', @(h,e)keyPressFcn_thresh_(hFig, 'n'));
% uimenu(mh_view,'Label', '.prm file', 'Callback', @edit_prm_);
uimenu(mh_view,'Label', 'Reset window positions', 'Callback', @reset_position_);

mh_proj = uimenu(hFig,'Label','Projection'); 
uimenu(mh_proj, 'Label', 'vpp', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'vpp', 'vmin'}));
uimenu(mh_proj, 'Label', 'pca', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'pca'}));
uimenu(mh_proj, 'Label', 'cov', 'Callback', @(h,e)proj_view_(h), ...
    'Checked', if_on_off_(P.vcFet_show, {'cov', 'spacetime'}));

mh_plot = uimenu(hFig,'Label','Plot'); 
uimenu(mh_plot, 'Label', 'All unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_);
uimenu(mh_plot, 'Label', 'Selected unit firing rate vs. aux. input', 'Callback', @(h,e)plot_aux_rate_(1));

mh_info = uimenu(hFig,'Label','','Tag', 'mh_info'); 
uimenu(mh_info, 'Label', 'Annotate unit', 'Callback', @unit_annotate_);
uimenu(mh_info, 'Label', 'single', 'Callback', @(h,e)unit_annotate_(h,e,'single'));
uimenu(mh_info, 'Label', 'multi', 'Callback', @(h,e)unit_annotate_(h,e,'multi'));
uimenu(mh_info, 'Label', 'noise', 'Callback', @(h,e)unit_annotate_(h,e,'noise'));
uimenu(mh_info, 'Label', 'clear annotation', 'Callback', @(h,e)unit_annotate_(h,e,''));
uimenu(mh_info, 'Label', 'equal to', 'Callback', @(h,e)unit_annotate_(h,e,'=%d'));

mh_history = uimenu(hFig,'Label', 'History', 'Tag', 'History'); 

mh_help = uimenu(hFig,'Label','Help'); 
uimenu(mh_help, 'Label', '[H]elp', 'Callback', @help_FigWav_);
uimenu(mh_help, 'Label', 'Wiki on GitHub', 'Callback', @(h,e)wiki_());
uimenu(mh_help, 'Label', 'About', 'Callback', @(h,e)msgbox_(about_()));
uimenu(mh_help, 'Label', 'Post an issue on GitHub', 'Callback', @(h,e)issue_('search'));
uimenu(mh_help, 'Label', 'Search issues on GitHub', 'Callback', @(h,e)issue_('post'));

drawnow;
set(hFig, 'OuterPosition', posvec);
end %func


%--------------------------------------------------------------------------
function S_fig = plot_spkwav_(S_fig, P)
global tnWav_spk tnWav_raw 
% fPlot_raw = 0;
if nargin<2, P = get0_('P'); end
[viSite_spk, S_clu] = get0_('viSite_spk', 'S_clu');
[cvrX, cvrY, cviSite] = deal(cell(S_clu.nClu, 1));
vnSpk = zeros(S_clu.nClu, 1);
miSites_clu = P.miSites(:, S_clu.viSite_clu);
if isfield(S_fig, 'maxAmp')
    maxAmp = S_fig.maxAmp;
else
    maxAmp = P.maxAmp;
end
for iClu = 1:S_clu.nClu        
    try
        viSpk_show = randomSelect_(S_clu_viSpk_(S_clu, iClu, viSite_spk), P.nSpk_show);
        if P.fWav_raw_show
            trWav1 = raw2uV_(tnWav_raw(:,:,viSpk_show), P);
            trWav1 = fft_lowpass_(trWav1, get_set_(P, 'fc_spkwar_show', []), P.sRateHz);
        else
            trWav1 = tnWav2uV_(tnWav_spk(:,:,viSpk_show), P);
        end        
        viSite_show = miSites_clu(:, iClu);
        [cvrY{iClu}, cvrX{iClu}] = tr2plot_(trWav1, iClu, viSite_show, maxAmp, P);
        cviSite{iClu} = viSite_show;
        vnSpk(iClu) = size(trWav1, 3); %subsample 
    catch
        disperr_();
    end
end
S = makeStruct_(cvrY, cviSite, vnSpk);
try
    set(S_fig.hSpkAll, 'XData', cell2mat_(cvrX), 'YData', cell2mat_(cvrY), 'UserData', S);
catch
    S_fig.hSpkAll = plot(S_fig.hAx, cell2mat_(cvrX), cell2mat_(cvrY), 'Color', [.5 .5 .5], 'LineWidth', .5); %, P.LineStyle); 
    set(S_fig.hSpkAll, 'UserData', S);
end
end %func


%--------------------------------------------------------------------------
function exit_manual_(src, event)
try    
    if ~ishandle(src), return; end
    if ~isvalid(src), return; end
    S0 = get(0, 'UserData'); 
    P = S0.P;
%     if ~get_set_([], 'fDebug_ui', 0)
    fExit = save_manual_(P);
    if ~fExit, return; end 
    if ~isfield(S0, 'csFig')
        S0.csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
    end
    delete_multi_(get_fig_all_(S0.csFig), src);
    close_(get_fig_('FigTrial'));
    close_(get_fig_('FigAux'));
catch
    disperr_();
    close(src);
end
end %func


%--------------------------------------------------------------------------
function fExit = save_manual_(varargin) 
% TODO: remove figure handles from S0
if nargin==1
    P = varargin{1};
else
    P = get0_('P');
end
vcFile_jrc = subsFileExt_(P.vcFile_prm, '_jrc.mat');
fExit = 1;
switch lower(questdlg_(['Save to ', vcFile_jrc, ' ?'], 'Confirmation', 'Yes'))
    case 'yes'
        hMsg = msgbox_('Saving');
%         assignWorkspace_(Sclu);
        save0_(vcFile_jrc);
        fExit = 1;
        close_(hMsg);
    case 'no'
        fExit = 1;
    case 'cancel' 
        fExit = 0;
        return;
end
end %func;


%--------------------------------------------------------------------------
function cvhHide_mouse = mouse_hide_(hFig, hObj_hide, S_fig)
% hide during mouse pan to speed up     
if nargin<3, S_fig = get(hFig, 'UserData'); end
% if nargin<3, S0 = get(0, 'UserData'); end
if nargin == 0 %clear field
%     try S_fig = rmfield(S_fig, 'vhFig_mouse'); catch; end
    try S_fig = rmfield(S_fig, 'cvhHide_mouse'); catch; end
else
    if ~isfield(S_fig, 'vhFig_mouse') && ~isfield(S_fig, 'cvhHide_mouse')
%         S_fig.vhFig_mouse = hFig;
        S_fig.cvhHide_mouse = {hObj_hide};    
    else
%         S_fig.vhFig_mouse(end+1) = hFig;
        S_fig.cvhHide_mouse{end+1} = hObj_hide;
    end
end
cvhHide_mouse = S_fig.cvhHide_mouse;
if nargout==0, set(hFig, 'UserData', S_fig); end
end %func


%--------------------------------------------------------------------------
function S_fig = plot_tnWav_clu_(S_fig, P)
% Substituting plot_spk_
S0 = get(0, 'UserData'); 
S_clu = S0.S_clu;
if ~isfield(P, 'LineWidth'), P.LineWidth=1; end
trWav_clu = ifeq_(P.fWav_raw_show, S_clu.tmrWav_raw_clu, S_clu.tmrWav_clu);
[nSamples, nSites, nClu] = size(trWav_clu);
nChans_show = size(P.miSites, 1);
miSites_clu = P.miSites(:, S_clu.viSite_clu);
% nSites = numel(P.viSite2Chan);

% determine x
x_offset = P.spkLim(2) / (diff(P.spkLim)+1); %same for raw and filt
vrX = (1:nSamples*nClu)/nSamples + x_offset;
vrX(1:nSamples:end) = nan;
vrX(nSamples:nSamples:end) = nan;
trWav_clu = trWav_clu / S_fig.maxAmp;

% nChans_show = size(P.miSites,1);
mrX = repmat(vrX(:), [1, nChans_show]);
mrX = reshape(mrX, [nSamples, nClu, nChans_show]);
mrX = reshape(permute(mrX, [1 3 2]), [nSamples*nChans_show, nClu]);

mrY = zeros(nSamples * nChans_show, nClu, 'single');
for iClu=1:nClu
    viSites1 = miSites_clu(:,iClu);
    mrY1 = trWav_clu(:,viSites1,iClu);
    mrY1 = bsxfun(@plus, mrY1, single(viSites1'));
    mrY(:,iClu) = mrY1(:);
end
    
% if isempty(P.LineStyle)
if isfield(S_fig, 'vhPlot')
    plot_update_(S_fig.vhPlot, mrX, mrY); 
else
    S_fig.vhPlot = plot_group_(S_fig.hAx, mrX, mrY, 'LineWidth', P.LineWidth); 
end
% else
%     S_fig.vhPlot = plot_group_(S_fig.hAx, mrX, mrY, P.LineStyle, 'LineWidth', P.LineWidth); 
% end
set(S_fig.hAx, 'YTick', 1:nSites, 'XTick', 1:nClu);
end %func


%--------------------------------------------------------------------------
function [hFig, S_fig] = plot_FigWavCor_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
S_clu = S0.S_clu; P = S0.P;
[hFig, S_fig] = get_fig_cache_('FigWavCor'); 

figure_wait_(1, hFig);
nClu = S_clu.nClu;
% Plot
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.1 .1 .8 .8], 'XLimMode', 'manual', 'YLimMode', 'manual', 'Layer', 'top');
    set(S_fig.hAx, {'XTick', 'YTick'}, {1:nClu, 1:nClu});
    axis(S_fig.hAx, [0 nClu 0 nClu]+.5);
    axis(S_fig.hAx, 'xy');
    grid(S_fig.hAx, 'on');
    xlabel(S_fig.hAx, 'Clu#'); 
    ylabel(S_fig.hAx, 'Clu#');
    S_fig.hImWavCor = imagesc(S_clu.mrWavCor, P.corrLim);  %clears title and current figure
    S_fig.hCursorV = line([1 1], [.5 nClu+.5], 'Color', [0 0 0], 'LineWidth', 1.5); 
    S_fig.hCursorH = line([.5 nClu+.5], [1 1], 'Color', [1 0 0], 'LineWidth', 1.5);             
    colorbar(S_fig.hAx);
    S_fig.vcTitle = '[S]plit; [M]erge; [D]elete';
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigWavCor_);
    mouse_figure(hFig, S_fig.hAx, @button_FigWavCor_);
    S_fig.hDiag = plotDiag_([0, nClu, .5], 'Color', [0 0 0], 'LineWidth', 1.5);
else
    set(S_fig.hImWavCor, 'CData', S_clu.mrWavCor);
    set(S_fig.hCursorV, 'xdata', [1 1], 'ydata', [.5 nClu+.5]);
    set(S_fig.hCursorH, 'xdata', .5+[0 nClu], 'ydata', [1 1]);
end
% output
set(hFig, 'UserData', S_fig);
figure_wait_(0, hFig);
end %func


%--------------------------------------------------------------------------
% function mrWavCor = calc_WavCor_(S_clu, P)
% % calculate waveform correlation
% % update waveforms
% global tnWav_spk tnWav_raw
% S0 = get0_();
% 
% S_clu = clu2wav_(S_clu, S0.viSite_spk, tnWav_spk, tnWav_raw);
% mrWavCor = S_clu_wavcor_(S_clu, P);
% nClu = double(S_clu.nClu);
% mrWavCor(sub2ind(nClu * [1,1], 1:nClu, 1:nClu)) = clu_self_corr_(S0, tnWav_raw);
% set0_(S_clu, mrWavCor);
% % update selfcorr_clu
% end %func


%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr_(S_clu, tnWav_spk, iClu1)
% plot top half vs bottom half correlation. sum of vpp

% S_clu = S0.S_clu;
viSite_spk = get0_('viSite_spk');
if nargin < 3
    fprintf('Computing self correlation\n\t'); t1=tic;
    selfcorr = zeros(1, S_clu.nClu);
    for iClu=1:S_clu.nClu
        selfcorr(iClu) = S_clu_self_corr__(S_clu, tnWav_spk, iClu, viSite_spk);
        fprintf('.');
    end
    fprintf('\n\ttook %0.1fs\n', toc(t1));
else
    selfcorr = S_clu_self_corr__(S_clu, tnWav_spk, iClu1, viSite_spk);
end
end %func


%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr__(S_clu, tnWav_spk, iClu1, viSite_spk)
% cluster self-correlation. low means bad. return 1-corr score
MAX_SAMPLE = 4000;
if nargin<4, viSite_spk = get0_('viSite_spk'); end

[viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, viSite_spk);

viSpk_clu1 = randomSelect_(viSpk_clu1, MAX_SAMPLE);
% trWav1 = meanSubt_(single(tnWav_spk(:,:,viSpk_clu1))); 
trWav1 = tnWav_spk(:,:,viSpk_clu1);
vrVpp = squeeze_(squeeze_(max(trWav1(:,1,:)) - min(trWav1(:,1,:))));
% vrVpp = sum(squeeze_(max(tnWav_spk) - min(trWav1)));
[~, viSrt] = sort(vrVpp);
imid = round(numel(viSrt)/2);
mrWavA = meanSubt_(mean(trWav1(:, :, viSrt(1:imid)), 3));
mrWavB = meanSubt_(mean(trWav1(:, :, viSrt(imid+1:end)), 3));
% selfcorr = calc_corr_(mrWavA(:), mrWavB(:));
% selfcorr = mean(mean(zscore_(mrWavA) .* zscore_(mrWavB)));
% selfcorr = mean(zscore_(mrWavA(:)) .* zscore_(mrWavB(:)));
selfcorr = corr_(mrWavA(:), mrWavB(:));
end %func


%--------------------------------------------------------------------------
function [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S_clu, iClu1, viSite_spk)
% get a subset of cluster that is centered
% return only centered spikes
% if nargin<2, S0 = get(0, 'UserData'); end
% S_clu = S0.S_clu;
if nargin<3, viSite_spk = get0_('viSite_spk'); end
iSite_clu1 = S_clu.viSite_clu(iClu1);
viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
viSite_clu1 = viSite_spk(viSpk_clu1);
viiSpk_clu1 = find(viSite_clu1 == iSite_clu1);
viSpk_clu1 = viSpk_clu1(viiSpk_clu1);
end %func


%--------------------------------------------------------------------------
function vrX = wav_clu_x_(iClu, P)
% determine x range of a cluster
if P.fWav_raw_show
    nSamples = diff(P.spkLim*2) + 1;
    x_offset = (P.spkLim(2)*2) / nSamples + iClu - 1;
    vrX = (0:nSamples-1) / nSamples + x_offset;
else
    nSamples = diff(P.spkLim) + 1;
    x_offset = (P.spkLim(2)) / nSamples + iClu - 1;
    vrX = (1:nSamples) / nSamples + x_offset;
end
vrX([1,end]) = nan;
end %func


%--------------------------------------------------------------------------
function update_plot_(hPlot, vrX, vrY, S_plot)
% update the plot with new x and y 

if nargin<4, S_plot = []; end
if isempty(hPlot), return; end
% selective plot to speed up plotting speed
if isempty(vrY) || isempty(vrX)
    hide_plot_(hPlot);
%     set(hPlot, 'XData', nan, 'YData', nan);  %visible off
    return;
end

% only update if both x and y are changed
vrX1 = get(hPlot, 'XData');
vrY1 = get(hPlot, 'YData');
fUpdate = 1;
if (numel(vrX1) == numel(vrX)) && (numel(vrY1) == numel(vrY))
    if (std(vrX1(:) - vrX(:)) == 0) && (std(vrY1(:) - vrY(:)) == 0)
        fUpdate = 0; 
    end
end
if fUpdate, set(hPlot, 'xdata', vrX, 'ydata', vrY); end
if ~isempty(S_plot), set(hPlot, 'UserData', S_plot); end
end %func


%--------------------------------------------------------------------------
% find intersection of two limit ranges
function xlim1 = trim_lim_(xlim1, xlim0)
dx = diff(xlim1);

if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
xlim1(1) = max(xlim1(1), xlim0(1));
xlim1(2) = min(xlim1(2), xlim0(2));
end %func


%--------------------------------------------------------------------------
% 8/9/17 JJJ: Generalized to any figure objects
function S0 = keyPressFcn_cell_(hObject, csKey, S0)
% Simulate key press function

if nargin<3, S0 = get(0, 'UserData'); end
figure_wait_(1); 
event1.Key = '';
if ischar(csKey), csKey = {csKey}; end
nKeys = numel(csKey);
keyPressFcn_ = get(hObject, 'KeyPressFcn');
for i=1:nKeys
    event1.Key = csKey{i};
    S0 = keyPressFcn_(hObject, event1, S0);
end
% drawnow;
figure_wait_(0);
if nargout==0, set(0, 'UserData', S0); end
end %func


%--------------------------------------------------------------------------
function S0 = keyPressFcn_FigWav_(hObject, event, S0) %amp dist
% global mrWav
global fDebug_ui

if nargin<3, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
P.LineStyle=[];
nSites = numel(P.viSite2Chan);
hFig = hObject;
S_fig = get(hFig, 'UserData');

switch lower(event.Key)
    case {'uparrow', 'downarrow'}
        rescale_fig_(event, S0, P);
        clu_info_(S0); %update figpos

    case {'leftarrow', 'rightarrow'}
        % switch the current clu
        if ~key_modifier_(event, 'shift');
            if strcmpi(event.Key, 'leftarrow')
                if S0.iCluCopy == 1, return; end
                S0.iCluCopy = S0.iCluCopy - 1;
            else
                if S0.iCluCopy == S_clu.nClu, return; end
                S0.iCluCopy = S0.iCluCopy + 1;
            end
%             S0.iCluPaste = []; 
        else
            if isempty(S0.iCluPaste)
                S0.iCluPaste = S0.iCluCopy;
            end
            if strcmpi(event.Key, 'leftarrow')
                if S0.iCluPaste == 1, return; end
                S0.iCluPaste = S0.iCluPaste - 1;
            else
                if S0.iCluPaste == S_clu.nClu, return; end
                S0.iCluPaste = S0.iCluPaste + 1;
            end
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0); %select first clu
        
    case 'm', S0 = ui_merge_(); % merge clusters
        
%     case 'l' %LFP
%         error('LFP is not implemented');
%         vcFile_lfp = subsFileExt(P.vcFile, '.lfp');
%         mrLfp = loadWavCache(vcFile_lfp, 'mrLfp', P);
%         mrLfp = bsxfun(@minus, mrLfp, median(mrLfp, 2));
%         viTime1 = Sclu_viTime_(S0.S_clu, S0.iCluCopy);
%         lim_lfp = [-200, 200];
%         tr1 = single(mr2tr2(mrLfp, viTime1, lim_lfp));
%         mu1 = mean(tr1, 3);
%         sd1 = std(tr1, 0, 3);
%         vrT1 = (lim_lfp(1):lim_lfp(2)) / (P.sRateHz / P.nSkip_lfp / 1000);
%         iSite1 = S0.S_clu.viSite_clu(S0.iCluCopy);
%         figure(2001); errorbar(vrT1, mu1(:,iSite1), sd1(:,iSite1)); xlabel('Time (ms)');
%         title_(sprintf('Clu %d, n=%d', S0.iCluCopy, numel(viTime1)));
%         S0.S_clu. S0.iCluCopy
        
    case 'space'
        % auto-select nearest cluster for black
        mrWavCor = S_clu.mrWavCor;
        mrWavCor(S0.iCluCopy,S0.iCluCopy) = -inf;
        [~,S0.iCluPaste] = max(mrWavCor(:,S0.iCluCopy));
        set(0, 'UserData', S0);
        button_CluWav_simulate_([], S0.iCluPaste);
        
    case 's', auto_split_(1);
        
    case 'r' %reset view
        figure_wait_(1);
        axis([0, S0.S_clu.nClu + 1, 0, numel(P.viSite2Chan) + 1]);
        figure_wait_(0);
        
    case {'d', 'backspace', 'delete'}, S0 = ui_delete_();
        
    case 'z' %zoom
        iClu = S0.iCluCopy;
        iSiteClu = S_clu.viSite_clu(S0.iCluCopy);
        set_axis_(hFig, iClu+[-1,1]*6, iSiteClu+[-1,1]*(P.maxSite*2+1), [0 S_clu.nClu+1], [0 nSites+1]);

    case 'c', plot_FigCorr_(S0);
        
    case 'o' %overlap waveforms across sites
        hFig_temp = figure; hold on;
        if isempty(S0.iCluPaste)
            viSite1 = P.miSites(:,S_clu.viSite_clu(S0.iCluCopy));
            plot(zscore_(S_clu.tmrWav_clu(:,viSite1,S0.iCluCopy)));
        else
            viSite1 = P.miSites(:,S_clu.viSite_clu(S0.iCluCopy));            
            plot(zscore_(S_clu.tmrWav_clu(:,viSite1,S0.iCluCopy)), 'k');
            plot(zscore_(S_clu.tmrWav_clu(:,viSite1,S0.iCluPaste)), 'r');
        end
        title_(sprintf('Clu%d', S0.iCluCopy));
        grid on; 
%         if get_set_([], 'fDebug_ui', 0), close_(hFig_temp); end
        if fDebug_ui==1, close_(hFig_temp); end
        
    case 'v', plot_FigIsi_(S0);
        
    case 'a', update_spikes_(S0); clu_info_(S0);
        
    case 'f', clu_info_(S0);       
        
    case 'h', msgbox_(S_fig.csHelp, 1);

    case 'w', toggleVisible_(S_fig.hSpkAll); %toggle spike waveforms
        
    case 't', plot_FigTime_(S0); % time view

%     case 'q' %quality score (klustaviwa)
% %         https://github.com/klusta-team/klustaviewa/blob/ff7d1e4dcf4e2320e9e847382bb1d6cd2c083212/klustaviewa/stats/quality.py
%         if isempty(S0.hCopy), return ;end
%         set(hObject, 'Pointer', 'watch');
%         viTime1 = Sclu.cviTime_clu{S0.iCluCopy}; %(Sclu.viClu == S.iCluCopy);        
%         tmrWav1 = getSpkWav_(S0.mrWav, viTime1, P);
%         mrWav1_mu = mean(tmrWav1,3);
%         [~, iChMin] = min(mrWav1_mu(1-P.spkLim(1),:));
%         viChan1 = Sclu.P.miSites(:,iChMin);
%         mrFet1 = squeeze_(tmrWav1(1-P.spkLim(1),viChan1,:)); %sub-select chan
% %         mrFet1 = reshape(tmrWav1, [], size(tmrWav1,3));
%         mrDot1 = mrFet1' * mrFet1;
%         mrDot1 = triu(mrDot1);
%         vrDist1 = mrDot1(mrDot1~=0);
%         [nsamples, nchannels, nspikes] = size(tmrWav1(:,viChan1,:));
%         q = 1. / nsamples * max(squeeze_(mean(sum(tmr_ .^ 2), 3)));
%         fprintf('clu%d quality: %f\n', S0.iCluCopy, q);
        
    case 'j', plot_FigProj_(S0); %projection view
        
    case 'n'        
        if strcmpi(get(S_fig.vhText(1), 'Visible'), 'on')
            set(S_fig.vhText, 'Visible', 'off');
        else
            set(S_fig.vhText, 'Visible', 'on');
        end       
        
    case 'i', plot_FigHist_(S0); %ISI histogram
        
%     case 'f' %fano factor
%         if isempty(S.hCopy), return ;end        
%         vrTime1 = double(Sclu.viTime(Sclu.viClu==S.iCluCopy)) / Sclu.P.sRateHz;
%         fanoFactor(vrTime1);
%         try close(hMsgbox); catch, end
        
    case 'e', plot_FigMap_(S0);
        
    case 'u', update_FigCor_(S0);
        
    case 'p' %PSTH plot
        if isempty(P.vcFile_trial), msgbox_('''vcFile_trial'' not set.'); return; end
        plot_raster_(S0.P, [S0.iCluCopy, S0.iCluPaste], S0.S_clu); %psth
        
    otherwise, figure_wait_(0); %stop waiting
end

figure(hObject); %change the focus back to the current object

end %func


%--------------------------------------------------------------------------
function export_rate_()
S_clu = get0_('S_clu');
mrRate_clu = clu_rate_(S_clu); 
csMsg = assignWorkspace_(mrRate_clu);
fprintf(csMsg);
msgbox_(csMsg, 1);
end %func


%--------------------------------------------------------------------------
function plot_FigCorr_(S0)
% hFigCorr plot
jitter_ms = .5; % bin size for correlation plot
nLags_ms = 25; %show 25 msec

if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
P.jitter_ms = jitter_ms;
P.nLags_ms = nLags_ms;

[hFig, S_fig] = get_fig_cache_('FigCorr'); 
iClu1 = S0.iCluCopy;
iClu2 = S0.iCluPaste;
if isempty(iClu2), iClu2 = iClu1; end

jitter = round(P.sRateHz / 1000 * P.jitter_ms); %0.5 ms
nLags = round(P.nLags_ms / P.jitter_ms);

vi1 = int32(double(S_clu_time_(S_clu, iClu1)) /jitter);

if iClu1~=iClu2
    vi1 = [vi1, vi1-1, vi1+1]; %allow missing one
end
vi2 = int32(double(S_clu_time_(S_clu, iClu2)) /jitter);
viLag = -nLags:nLags;
vnCnt = zeros(size(viLag));
for iLag=1:numel(viLag)
    if iClu1 == iClu2 && viLag(iLag)==0, continue; end
    vnCnt(iLag) = numel(intersect(vi1, vi2+viLag(iLag)));
end
vrTime_lag = viLag * P.jitter_ms;

%--------------
% draw
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    S_fig.hBar = bar(vrTime_lag, vnCnt, 1);     
    xlabel('Time (ms)'); 
    ylabel('Counts');
    grid on; 
    set(S_fig.hAx, 'YScale', 'log');
else
    set(S_fig.hBar, 'XData', vrTime_lag, 'YData', vnCnt);
end
title_(S_fig.hAx, sprintf('Clu%d vs Clu%d', iClu1, iClu2));
xlim(S_fig.hAx, [-nLags, nLags] * P.jitter_ms);
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function plot_FigTime_(S0)
% plot FigTime window. Uses subsampled data

if nargin<1, S0 = get(0, 'UserData'); end
S_clu = S0.S_clu; P = S0.P; 
[hFig, S_fig] = get_fig_cache_('FigTime');

%----------------
% collect info
iSite = S_clu.viSite_clu(S0.iCluCopy);
[vrFet0, vrTime0] = getFet_site_(iSite);    % plot background    
[vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy); % plot iCluCopy

vcTitle = '[H]elp; (Sft)[Left/Right]:Sites/Features; (Sft)[Up/Down]:Scale; [B]ackground; [S]plit; [R]eset view; [P]roject; [M]erge; (sft)[Z] pos; [E]xport selected; [C]hannel PCA';
if ~isempty(S0.iCluPaste)
    [vrFet2, vrTime2] = getFet_site_(iSite, S0.iCluPaste);
    vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', S0.iCluCopy, S0.iCluPaste, vcTitle);
else
    vrFet2 = [];
    vrTime2 = [];
    vcTitle = sprintf('Clu%d (black); %s', S0.iCluCopy, vcTitle);
end
time_lim = double([0, S0.viTime_spk(end)] / P.sRateHz);

%------------
% draw
if isempty(S_fig)
    S_fig.maxAmp = P.maxAmp;
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');
    
    % first time
    S_fig.hPlot0 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(1,:), 'MarkerSize', 5, 'LineStyle', 'none');
    S_fig.hPlot1 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(2,:), 'MarkerSize', 5, 'LineStyle', 'none');
    S_fig.hPlot2 = line(nan, nan, 'Marker', '.', 'Color', P.mrColor_proj(3,:), 'MarkerSize', 5, 'LineStyle', 'none');   %place holder  
    xlabel('Time (s)'); 
    grid on;    
    
    % rectangle plot
    vrPos_rect = [time_lim(1), S_fig.maxAmp, diff(time_lim), S_fig.maxAmp];
    S_fig.hRect = imrect_(S_fig.hAx, vrPos_rect); %default position?
    if ~isempty(S_fig.hRect)
        setColor(S_fig.hRect, 'r');
        setPositionConstraintFcn(S_fig.hRect, ...
            makeConstrainToRectFcn('imrect',time_lim, [-4000 4000]));  
    end
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigTime_);
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
    if ~isempty(P.time_tick_show) %tick mark
        set(S_fig.hAx, 'XTick', time_lim(1):P.time_tick_show:time_lim(end));
    end
%     S_fig.iFet = 1;
%     S_fig.csFet = {'vpp'};
%     for iFet = 1:P.nPcPerChan
%         S_fig.csFet{iFet+1} = sprintf('%s%d', P.vcFet, iFet);
%     end
%     S_fig.csFet = S_fig.csFet(1:1+P.nPcPerChan);
end
vpp_lim = [0, S_fig.maxAmp];
% iFet = S_fig.iFet;
% iFet = 1;
if ~isfield(S_fig, 'iSite'), S_fig.iSite = []; end
update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
update_plot_(S_fig.hPlot1, vrTime1, vrFet1);
update_plot_(S_fig.hPlot2, vrTime2, vrFet2);
imrect_set_(S_fig.hRect, time_lim, vpp_lim);
mouse_figure(hFig, S_fig.hAx); % allow zoom using wheel
% button click function to select individual spikes, all spikes plotted

if isfield(S_fig, 'vhAx_track')
    toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 0);
    toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot1, S_fig.hPlot2, S_fig.hPlot0}, 1);
end

if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);

axis(S_fig.hAx, [time_lim, vpp_lim]);
title_(S_fig.hAx, vcTitle);    
ylabel(S_fig.hAx, vcYlabel);

S_fig = struct_merge_(S_fig, makeStruct_(iSite, time_lim, P, vpp_lim, viSpk1));
S_fig.csHelp = {...
    'Up/Down: change channel', ...
    'Left/Right: Change sites', ...
    'Shift + Left/Right: Show different features', ...
    'r: reset scale', ...
    'a: auto-scale', ...
    'c: show pca across sites', ...
    'e: export cluster info', ...
    'f: export cluster feature', ...
    'Zoom: mouse wheel', ...
    'H-Zoom: press x and wheel. space to reset', ...
    'V-Zoom: press y and wheel. space to reset', ...
    'Drag while pressing wheel: pan'};
        
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, iClu)
% just specify iSite to obtain background info
% 2016 07 07 JJJ
% return feature correspojnding to a site and cluster
% requiring subsampled info: cvrVpp_site and cmrFet_site. store in S0
% global mrWav

S0 = get(0, 'UserData');
% S_clu = S0.S_clu;
P = S0.P;
if ~isfield(P, 'vcFet_show'), P.vcFet_show = 'vpp'; end
% vrFet1 = []; vrTime1=[];
% iFet = str2double(P.vcFet_show(end));
if nargin < 2, iClu = []; end
[vrFet1, viSpk1] = getFet_clu_(iClu, iSite);
vrTime1 = double(S0.viTime_spk(viSpk1)) / P.sRateHz;
        
% if nargin < 2
%    % get background feature      
%    switch lower(P.vcFet_show)
%        case {'vmin', 'vpp'}
%             vrFet1 = S0.cvrVpp_site{iSite};   
%        otherwise
%            error('not implemented yet'); % list of background spikes
% %                 vrFet1 = S0.cmrFet_site{iSite}(:,iFet);
%     end
%    vrTime1 = S0.cvrTime_site{iSite};
%    viSpk1 = []; %@TODO: Sclu.cviSpk_site{iSite};
% else        
%     [vrFet1, viSpk1] = getFet_clu_(iClu, iSite);
%     vrTime1 = double(S0.viTime_spk(viSpk1)) / P.sRateHz;
% end

% label
switch lower(P.vcFet_show)
    case {'vpp', 'vmin'} %voltage feature
        vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.vcFet_show);
    otherwise %other feature options
        vcYlabel = sprintf('Site %d (%s)', iSite, P.vcFet_show);
end

end %func 


%--------------------------------------------------------------------------
function [viTime1, viSpk1, viSpk2] = S_clu_time_(S_clu, iClu)
% return time of cluster time in adc sample index unit.
% viSpk2: spike indices of centered spikes in cluster iClu

S0 = get(0, 'UserData');
if isfield(S_clu, 'cviSpk_clu')
    viSpk1 = S_clu.cviSpk_clu{iClu};
    viTime1 = S0.viTime_spk(viSpk1);
else
    viSpk1 = find(S_clu.viClu == iClu);
    viTime1 = S0.viTime_spk(viSpk1);
end
if nargout>=3
    iSite1 = S_clu.viSite_clu(iClu);
    viSpk2 = viSpk1(S0.viSite_spk(viSpk1) == iSite1);
end
end %func


%--------------------------------------------------------------------------
function imrect_set_(hRect, xpos, ypos)
vrPos = getPosition(hRect);
if ~isempty(xpos)
    vrPos(1) = min(xpos);
    vrPos(3) = abs(diff(xpos));
end
if ~isempty(ypos)
    vrPos(2) = min(ypos);
    vrPos(4) = abs(diff(ypos));
end
setPosition(hRect, vrPos);
end %func


%--------------------------------------------------------------------------
function flag = isVisible_(hObj)
flag = strcmpi(get(hObj, 'Visible'), 'on');
end %func


%--------------------------------------------------------------------------
function plot_FigProj_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
S_clu = S0.S_clu; P = S0.P;
[hFig, S_fig] = get_fig_cache_('FigProj');

iClu1 = S0.iCluCopy;
iClu2 = S0.iCluPaste;
update_plot2_proj_(); %erase prev objects

%---------------
% Compute
iSite1 = S_clu.viSite_clu(iClu1);
% miSites = P.miSites;
if ~isfield(P, 'viSites_show')
    P.viSites_show = sort(P.miSites(:, iSite1), 'ascend');
end
viSites_show = P.viSites_show;
nSites = numel(P.viSites_show);
cell_plot = {'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'};
switch lower(P.vcFet_show)
    case {'vpp', 'vmin', 'vmax'}
        vcXLabel = 'Site # (%0.0f \\muV; upper: V_{min}; lower: V_{max})';
        vcYLabel = 'Site # (%0.0f \\muV_{min})';
    otherwise
        vcXLabel = sprintf('Site # (%%0.0f %s; upper: %s1; lower: %s2)', P.vcFet_show, P.vcFet_show, P.vcFet_show);
        vcYLabel = sprintf('Site # (%%0.0f %s)', P.vcFet_show);    
end
vcTitle = '[H]elp; [S]plit; [B]ackground; (Sft)[Up/Down]:Scale; [Left/Right]:Sites; [M]erge; [F]eature';

%----------------
% display
if isempty(S_fig)
    S_fig.maxAmp = P.maxAmp;    
    S_fig.hAx = axes_new_(hFig);
    set(S_fig.hAx, 'Position', [.1 .1 .85 .85], 'XLimMode', 'manual', 'YLimMode', 'manual');
    S_fig.hPlot0 = line(nan, nan, 'Color', P.mrColor_proj(1,:), 'Parent', S_fig.hAx);
    S_fig.hPlot1 = line(nan, nan, 'Color', P.mrColor_proj(2,:), 'Parent', S_fig.hAx); %place holder
    S_fig.hPlot2 = line(nan, nan, 'Color', P.mrColor_proj(3,:), 'Parent', S_fig.hAx); %place holder
    set([S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2], cell_plot{:}); %common style
    S_fig.viSites_show = []; %so that it can update
    S_fig.vcFet_show = 'vpp';
    % plot boundary
    plotTable_([0, nSites], '-', 'Color', [.5 .5 .5]); %plot in one scoop
    plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5); %plot in one scoop
    mouse_figure(hFig);
    set(hFig, 'KeyPressFcn', @keyPressFcn_FigProj_);
    S_fig.cvhHide_mouse = mouse_hide_(hFig, S_fig.hPlot0, S_fig);
    set_fig_(hFig, S_fig);
end
% get features for x0,y0,S_plot0 in one go
%[mrMin, mrMax, vi0, vi1, vi2] = fet2proj_(S0, P.viSites_show);
[mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, P.viSites_show);
% S_fig.maxAmp %debug
if ~isfield(S_fig, 'viSites_show'), S_fig.viSites_show = []; end
if ~equal_vr_(S_fig.viSites_show, P.viSites_show) || ...
    ~equal_vr_(S_fig.vcFet_show, P.viSites_show)
    plot_proj_(S_fig.hPlot0, mrMin0, mrMax0, P, S_fig.maxAmp);
end

plot_proj_(S_fig.hPlot1, mrMin1, mrMax1, P, S_fig.maxAmp);
if ~isempty(iClu2)
    plot_proj_(S_fig.hPlot2, mrMin2, mrMax2, P, S_fig.maxAmp);
    vcTitle = sprintf('Clu%d (black), Clu%d (red); %s', iClu1, iClu2, vcTitle);
else
    update_plot_(S_fig.hPlot2, nan, nan);
    vcTitle = sprintf('Clu%d (black); %s', iClu1, vcTitle);
end

% Annotate axes
axis(S_fig.hAx, [0 nSites 0 nSites]);
set(S_fig.hAx,'XTick',.5:1:nSites,'YTick',.5:1:nSites, 'XTickLabel', P.viSites_show, 'YTickLabel', P.viSites_show, 'Box', 'off');
xlabel(S_fig.hAx, sprintf(vcXLabel, S_fig.maxAmp));   
ylabel(S_fig.hAx, sprintf(vcYLabel, S_fig.maxAmp));  
title_(S_fig.hAx, vcTitle);
vcFet_show = P.vcFet_show;
S_fig = struct_merge_(S_fig, ...
    makeStruct_(vcTitle, iClu1, iClu2, viSites_show, vcXLabel, vcYLabel, vcFet_show));
S_fig.csHelp = { ...
    '[D]raw polygon', ...
    '[S]plit cluster', ...
    '(shift)+Up/Down: change scale', ...
    '[R]eset scale', ...
    'Zoom: mouse wheel', ...
    'Drag while pressing wheel: pan'};
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function update_plot2_proj_(vrX, vrY)
if nargin==0, vrX=nan; vrY=nan; end
[hFig, S_fig] = get_fig_cache_('FigProj');
% erase polygon
if nargin==0
    try
        update_plot_(S_fig.hPlot2, vrX, vrY);
        delete(findobj(get(S_fig.hAx, 'Child'), 'Type', 'hggroup'));
    catch
        ;
    end
end
end


%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_refresh_(S_clu, fRemoveEmpty)

if nargin<2, fRemoveEmpty=1; end
nClu = double(max(S_clu.viClu));
S_clu.nClu = nClu;
viSite_spk = get0_('viSite_spk');
if isfield(S_clu, 'viSpk_shank')
    viSite_spk = viSite_spk(S_clu.viSpk_shank);
end
% gviClu = gpuArray_(S_clu.viClu);
% S_clu.cviSpk_clu = arrayfun(@(iClu)gather(find(gviClu==iClu)), 1:nClu, 'UniformOutput', 0);
S_clu.cviSpk_clu = arrayfun(@(iClu)find(S_clu.viClu==iClu), 1:nClu, 'UniformOutput', 0);
S_clu.vnSpk_clu = cellfun(@numel, S_clu.cviSpk_clu); 
S_clu.viSite_clu = double(arrayfun(@(iClu)mode(viSite_spk(S_clu.cviSpk_clu{iClu})), 1:nClu));
if fRemoveEmpty, [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu); end
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_map_index_(S_clu, viMap_clu)
% update viClu
vlPos = S_clu.viClu > 0;
viMap_clu = int32(viMap_clu);
S_clu.viClu(vlPos) = viMap_clu(S_clu.viClu(vlPos)); %translate cluster number
S_clu = S_clu_refresh_(S_clu, 0); % computational efficiency
% S_clu = S_clu_count_(S_clu);
end %func


%--------------------------------------------------------------------------
function [S_clu, vlKeep_clu] = S_clu_remove_empty_(S_clu)
% remove empty clusters, called by S_clu_refresh_ only

vlKeep_clu = S_clu.vnSpk_clu>0;
viClu_removed = find(~vlKeep_clu);
if isempty(viClu_removed), return; end

% fprintf('%d null-clusters removed\n', numel(viClu_removed));
% S_clu.vnSpk_clu(viClu_removed) = [];
% S_clu.cviSpk_clu(viClu_removed) = [];
% S_clu.viSite_clu(viClu_removed) = [];
% if isfield(S_clu, 'csNote_clu'), S_clu.csNote_clu = S_clu.csNote_clu(vlKeep_clu); end
if isfield(S_clu, 'mrWavCor'), S_clu.mrWavCor = S_clu.mrWavCor(vlKeep_clu, vlKeep_clu); end
% waveform
S_clu = S_clu_select_(S_clu, vlKeep_clu);
if min(S_clu.viClu) < 1
    S_clu.viClu(S_clu.viClu<1) = 0;
    [~,~,S_clu.viClu] = unique(S_clu.viClu+1);        
    S_clu.viClu = S_clu.viClu-1;
else
    [~,~,S_clu.viClu] = unique(S_clu.viClu);        
end
S_clu.viClu = int32(S_clu.viClu);
S_clu.nClu = double(max(S_clu.viClu));
end %func


%--------------------------------------------------------------------------
function [mrFet1, mrFet2, mrFet3, trWav2_spk] = trWav2fet_(tnWav1_spk, P, nSites_spk, viSite2_spk)
% [mrFet1, mrFet2, mrFet3, trWav_spk2] = trWav2fet_(tnWav_spk1, P)
% mrFet = trWav2fet_(tnWav_spk1, P)
if nargin<3, nSites_spk = []; end
if nargin<4, viSite2_spk = []; end
% global mrPv_global; % tnWav_spk

% if nargin<3, viSites_ref = spkwav_car_init_(P); end
[mrFet1, mrFet2, mrFet3] = deal(single([]));

% if strcmpi(P.vcFet, 'xpca')
%     [mrFet1, mrFet2] = pca_tr_(trWav_spk1);
%     mrFet1 = abs(mrFet1);
%     mrFet2 = abs(mrFet2);
%     mrFet3 = [];
%     return;
% end
fMeanSubt = 1;
% MAX_SAMPLE = 10000;
% viSites_ref = spkwav_car_init_(P);

%tnWav_spk1 = gpuArray_(tnWav_spk1, P.fGpu);
trWav2_spk = single(permute(tnWav1_spk, [1,3,2]));    
trWav2_spk = spkwav_car_(trWav2_spk, P, nSites_spk, viSite2_spk);

% [mrFet1, mrFet2, mrFet3] = tr2Vpp(trWav_spk1, P);
switch lower(P.vcFet) %{'xcor', 'amp', 'slope', 'pca', 'energy', 'vpp', 'diff248', 'spacetime'}       
    case {'spacetime', 'cov', 'cov2'}
        [mrFet1, mrFet2] = trWav2fet_cov_(trWav2_spk, P);
%         nDelay = 3;
%         if fMeanSubt
%             gtrWav1 = meanSubt_(trWav_spk1); 
%         else
%             gtrWav1 = (trWav_spk1); 
%         end
%         mr1 = gtrWav1(:,:,1);
%         mr1 = bsxfun(@rdivide, mr1, sqrt(mean(mr1.^2))); %zscore fast        
%         mr2 = mr1([ones(1,nDelay),1:end-nDelay],:,1); %time shift
%         mrFet1 = mean(gtrWav1 .* repmat(mr1, [1,1,size(gtrWav1,3)]), 1);
%         mrFet2 = mean(gtrWav1 .* repmat(mr2, [1,1,size(gtrWav1,3)]), 1);
%         mrFet1 = permute(mrFet1, [3,2,1]);
%         mrFet2 = permute(mrFet2, [3,2,1]);
%         if strcmpi(P.vcFet, 'cov2')
%             mrFet1(viSites_ref,:) = 0;
%             mrFet2(viSites_ref,:) = 0;
%         end
    case 'cov_prev'
        nDelay = 3;
        gtrWav1 = meanSubt_(trWav2_spk); 
        mr1 = zscore_(gtrWav1(:,:,1));                
        mr2 = zscore_(gtrWav1([ones(1,nDelay),1:end-nDelay],:,1)); 
        mrFet1 = mean(gtrWav1 .* repmat(mr1, [1,1,size(gtrWav1,3)]), 1);
        mrFet2 = mean(gtrWav1 .* repmat(mr2, [1,1,size(gtrWav1,3)]), 1);
        mrFet1 = shiftdim(mrFet1,1)';
        mrFet2 = shiftdim(mrFet2,1)';       
        
    case {'vpp', 'vppsqrt'}
        mrFet1 = shiftdim(max(trWav2_spk) - min(trWav2_spk))';
        if strcmpi(P.vcFet, 'vppsqrt'), mrFet1 = sqrt(mrFet1); end
        
    case {'amp', 'vmin'}
        mrFet1 = shiftdim(abs(min(trWav2_spk)))';
        
    case {'vminmax', 'minmax'}
        mrFet1 = shiftdim(abs(min(trWav2_spk)))';
        mrFet2 = shiftdim(abs(max(trWav2_spk)))';
        
    case 'energy'
        mrFet1 = shiftdim(std(trWav2_spk,1))';
        
    case 'energy2'
        nDelay = 3;
        mrFet1 = shiftdim(std(trWav2_spk,1))';
        trcov_ = @(a,b)shiftdim(sqrt(abs(mean(a.*b) - mean(a).*mean(b))));        
        mrFet2 = trcov_(trWav2_spk(1:end-nDelay,:,:), trWav2_spk(nDelay+1:end,:,:))';
        %mrFet1 = shiftdim(std(trWav_spk1,1))';
        
    case {'pca', 'gpca'}
        %  Compute PrinVec, 2D, max channel only
        if strcmpi(P.vcFet, 'pca')
            mrPv = tnWav2pv_(trWav2_spk, P);
        else
%             trWav_spk1 = spkwav_car_(trWav_spk1, viSites_ref);
            mrPv_global = get0_('mrPv_global');
            if isempty(mrPv_global)
                [mrPv_global, vrD_global] = tnWav2pv_(trWav2_spk, P);
                set0_(mrPv_global, vrD_global);
            end
            mrPv = mrPv_global;
        end
        
        % project
        dimm1 = size(trWav2_spk);
        if ismatrix(trWav2_spk), dimm1(end+1) = 1; end
%         mrWav_spk1 = meanSubt_(reshape(trWav2_spk, dimm1(1), []));
        mrWav_spk1 = reshape(trWav2_spk, dimm1(1), []);
        mrFet1 = reshape(mrPv(:,1)' * mrWav_spk1, dimm1(2:3))';
        if P.nPcPerChan >= 2
            mrFet2 = reshape(mrPv(:,2)' * mrWav_spk1, dimm1(2:3))';
        end
        if P.nPcPerChan >= 3
            mrFet3 = reshape(mrPv(:,3)' * mrWav_spk1, dimm1(2:3))';
        end
        
%     case 'std'
%         mrFet1 = squeeze_(std(trWav_spk1,1,2))';
end
if nargout==1
    mrFet1 = cat(1, mrFet1, mrFet2, mrFet3);
end
end %func


%--------------------------------------------------------------------------
function [mrPv, vrD1] = tnWav2pv_(tr, P)
%tr: nSamples x nSpikes x nChans
global tnWav_spk
MAX_SAMPLE = 10000;        

if nargin<2, P = get0_('P'); end
% if nargin<3, viSites_ref = []; end
if isempty(tr)
    nSpk = size(tnWav_spk,3);
    viSpk_sub = subsample_vr_(1:nSpk, MAX_SAMPLE);
    tr = permute(tnWav_spk(:,:,viSpk_sub), [1 3 2]);
    tr = single(tr);
    tr = spkwav_car_(tr, P);
    mrSpkWav1 = tr(:,:,1);
else
    viSpk_sub = subsample_vr_(1:size(tr,2), MAX_SAMPLE);
    mrSpkWav1 = tr(:,viSpk_sub, 1);    
end
% tr = single(tr);
% if ~isempty(viSites_ref), tr = spkwav_car_(tr, viSites_ref); end

% mrCov = meanSubt_(mrSpkWav1);
mrCov = mrSpkWav1 * mrSpkWav1';
[mrPv1, vrD1] = eig(mrCov);
mrPv1 = zscore_(fliplr(mrPv1)); % sort largest first
vrD1 = flipud(diag(vrD1));

% spike center should be negative
iMid = 1-P.spkLim(1);
vrSign = (mrPv1(iMid,:) < 0) * 2 - 1; %1 or -1 depending on the sign
mrPv = bsxfun(@times, mrPv1, vrSign);
end


%--------------------------------------------------------------------------
function [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, viSites0)
% show spikes excluding the clusters excluding clu1 and 2
P = S0.P;
S_clu = S0.S_clu;
iClu1 = S0.iCluCopy;
iClu2 = S0.iCluPaste;

% select subset of spikes
viSpk0 = find(ismember(S0.viSite_spk, viSites0));
viTime0 = S0.viTime_spk(viSpk0);
%time filter
if ~isfield(P, 'tlim_proj'), P.tlim_proj = []; end
if ~isempty(P.tlim_proj) 
    nlim_proj = round(P.tlim_proj * P.sRateHz);
    viSpk01 = find(viTime0>=nlim_proj(1) & viTime0<=nlim_proj(end));
    viSpk0 = viSpk0(viSpk01);
    viTime0 = viTime0(viSpk01);
end
viClu0 = S_clu.viClu(viSpk0);
viSpk00 = randomSelect_(viSpk0, P.nShow_proj*2);
viSpk01 = randomSelect_(viSpk0(viClu0 == iClu1), P.nShow_proj);
if ~isempty(iClu2)
    viSpk02 = randomSelect_(viSpk0(viClu0 == iClu2), P.nShow_proj);
else
    [mrMin2, mrMax2] = deal([]);
end
switch lower(P.vcFet_show)
    case {'pca'} %channel by channel pca. do it by channel
        % determine pca vector from cluster 1
        [mrPv1, mrPv2] = pca_pv_spk_(S_clu.cviSpk_clu{iClu1}, viSites0);
        [mrMin0, mrMax0] = pca_pc_spk_(viSpk00, viSites0, mrPv1, mrPv2); %getall spikes whose center lies in certain range
        [mrMin1, mrMax1] = pca_pc_spk_(viSpk01, viSites0, mrPv1, mrPv2); %getall spikes whose center lies in certain range
        if ~isempty(iClu2)  
            [mrMin2, mrMax2] = pca_pc_spk_(viSpk02, viSites0, mrPv1, mrPv2);
        end 
        
    otherwise % generic
        [mrMin0, mrMax0] = getFet_spk_(viSpk00, viSites0); %getall spikes whose center lies in certain range
        [mrMin1, mrMax1] = getFet_spk_(viSpk01, viSites0); %getall spikes whose center lies in certain range
        if ~isempty(iClu2)  
            [mrMin2, mrMax2] = getFet_spk_(viSpk02, viSites0);
        end            
end %switch
[mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = ...
    multifun_(@(x)abs(x), mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2);
end %func


%--------------------------------------------------------------------------
function [mrMin, mrMax] = getFet_spk_(viSpk1, viSites1)
% get feature for the spikes of interest
% global tnWav_spk

P = get0_('P');
switch lower(P.vcFet_show)
    case {'vmin', 'vpp'}   
        tnWav_spk1 = tnWav2uV_(tnWav_spk_sites_(viSpk1, viSites1));
        [mrMin, mrMax] = multifun_(@(x)abs(permute(x,[2,3,1])), min(tnWav_spk1), max(tnWav_spk1));
    case {'cov', 'spacetime'}
        [mrMin, mrMax] = calc_cov_spk_(viSpk1, viSites1);
    case 'pca'
        [mrMin, mrMax] = pca_pc_spk_(viSpk1, viSites1); %getall spikes whose center lies in certain range
    otherwise
        error('not implemented yet');
%         [mrFet1, mrFet2, mrFet3] = trWav2fet_(tnWav_spk1, P);
%         mrMin = (squeeze_(trFet_(1, viSites0, viSpk0)));
%         mrMax = (squeeze_(trFet_(2, viSites0, viSpk0)));
end
end %func


%--------------------------------------------------------------------------
function [mrVpp1, mrVpp2] = calc_cov_spk_(viSpk1, viSites1)
global tnWav_spk
% try
[viSite_spk, P] = get0_('viSite_spk', 'P');
% nT_spk = size(tnWav_spk, 1);
nSpk1 = numel(viSpk1);
viSites_spk1 = viSite_spk(viSpk1);
tnWav_spk1 = gpuArray_(tnWav_spk(:,:,viSpk1), P.fGpu); 
nSites_spk = 1 + P.maxSite * 2;
[mrVpp1_, mrVpp2_] = trWav2fet_(tnWav_spk1, P, nSites_spk);
[mrVpp1_, mrVpp2_] = multifun_(@(x)gather_(abs(x)), mrVpp1_, mrVpp2_);

% re-project to common basis
viSites_spk_unique = unique(viSites_spk1);
[mrVpp1, mrVpp2] = deal(zeros([numel(viSites1), nSpk1], 'like', mrVpp1_));
for iSite1 = 1:numel(viSites_spk_unique) %only care about the first site

        iSite11 = viSites_spk_unique(iSite1); %center sites group
        viSpk11 = find(viSites_spk1 == iSite11); %dangerous error
        viSites11 = P.miSites(:, iSite11);        
        [vlA11, viiB11] = ismember(viSites11, viSites1);
        mrVpp1(viiB11(vlA11),viSpk11) = mrVpp1_(vlA11,viSpk11);
        mrVpp2(viiB11(vlA11),viSpk11) = mrVpp2_(vlA11,viSpk11);
end    
% catch
%     disp('a');
% end
end %func


%--------------------------------------------------------------------------
function tnWav1 = tnWav1_sites_1_(tnWav1_, miSites1, viSites1)
[nT_spk, nSites_spk, nSpk1] = size(tnWav1_);
nSites1 = numel(viSites1);
% assert(nSites_spk==nSites1, 'tnWav1_sites_: nSites must agree');
tnWav1 = zeros([nT_spk, nSpk1, nSites1], 'like', tnWav1_); %subset of spk, complete
for iSite1 = 1:nSites1
    [viSite11, viiSpk11] = find(miSites1 == viSites1(iSite1));
    nSpk11 = numel(viiSpk11);
    mnWav_spk11 = reshape(tnWav1_(:, :, viiSpk11), nT_spk, []);
    mnWav_spk11 = mnWav_spk11(:, sub2ind([nSites_spk, nSpk11], viSite11', 1:nSpk11));
    tnWav1(:, viiSpk11, iSite1) = mnWav_spk11;
end
tnWav1 = permute(tnWav1, [1,3,2]);
end %func


%--------------------------------------------------------------------------
function tnWav_spk1 = tnWav_spk_sites_(viSpk1, viSites1, fWav_raw_show)
% reorder tnWav1 to viSites1
global tnWav_spk tnWav_raw
P = get0_('P');
% if nargin<3, fWav_raw_show = P.fWav_raw_show; end
if nargin<3, fWav_raw_show = 0; end
[viSite_spk, P] = get0_('viSite_spk', 'P');
if fWav_raw_show
    tnWav = tnWav_raw;
else
    tnWav = tnWav_spk;
end
nT_spk = size(tnWav, 1);
nSpk1 = numel(viSpk1);
viSites_spk1 = viSite_spk(viSpk1);
tnWav_spk1_ = gpuArray_(tnWav(:,:,viSpk1), P.fGpu);
viSites_spk_unique = unique(viSites_spk1);
tnWav_spk1 = zeros([nT_spk, numel(viSites1), nSpk1], 'like', tnWav_spk1_);
for iSite1 = 1:numel(viSites_spk_unique) %only care about the first site
    iSite11 = viSites_spk_unique(iSite1); %center sites group
    viSpk11 = find(viSites_spk1 == iSite11); %dangerous error
    viSites11 = P.miSites(:, iSite11);        
    [vlA11, viiB11] = ismember(viSites11, viSites1);
    tnWav_spk1(:,viiB11(vlA11),viSpk11) = tnWav_spk1_(:,vlA11,viSpk11);
end    
tnWav_spk1 = gather_(tnWav_spk1);
end %func


%--------------------------------------------------------------------------
function tnWav1 = tnWav1_sites_2_(tnWav1_, viSites_spk1, viSites1, P)
% reorder tnWav1 to viSites1

[nT_spk, nSites_spk, nSpk1] = size(tnWav1_);
% nSites1 = numel(viSites1);
nSites = numel(P.viSite2Chan);
tnWav1 = zeros([nT_spk, nSites, nSpk1], 'like', tnWav1_); %full
% assert(nSites_spk==nSites1, 'tnWav1_sites_: nSites must agree');
% tnWav1 = zeros([nT_spk, nSpk1, nSites1], 'like', tnWav1_); %subset of spk, complete
for iSite1 = 1:numel(viSites1)
    iSite11 = viSites1(iSite1);
    viSpk1 = find(viSites_spk1 == iSite11);
    if isempty(viSpk1), return; end
    viSites11 = P.miSites(:, iSite11);
    tnWav1(:,viSites11,viSpk1) = tnWav1_(:,:,viSpk1);
end
% tnWav1 = permute(tnWav1, [1,3,2]);
tnWav1 = tnWav1(:, viSites1, :); %reduced subset
end %func


%--------------------------------------------------------------------------
function flag = equal_vr_(vr1, vr2)
if all(size(vr1) == size(vr2))
    ml = vr1 == vr2;
    flag = all(ml(:));
else
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function plot_proj_(hPlot, mrMin, mrMax, P, maxAmp)
if nargin<5
    [hFig, S_fig] = get_fig_cache_('FigProj');
    maxAmp = S_fig.maxAmp;
end
[vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, P.maxSite_show);

% make struct
maxPair = P.maxSite_show;
viSites_show = P.viSites_show;
S_plot = makeStruct_(mrMax, mrMin, viSites_show, viPlot, tr_dim, maxPair, maxAmp);

update_plot_(hPlot, vrX, vrY, S_plot);
end %func


%--------------------------------------------------------------------------
function [vrX, vrY, viPlot, tr_dim] = amp2proj_(mrMin, mrMax, maxAmp, maxPair)
if nargin<4, maxPair = []; end
P = get0_('P');
% switch lower(P.vcFet_show)
%     case {'vpp', 'vmin', 'vmax'}
%         mrMax = linmap_(mrMax', [0, maxAmp/2], [0,1], 1);
%         mrMin = linmap_(mrMin', [0, maxAmp], [0,1], 1);
%     otherwise
        mrMax = linmap_(mrMax', [0, 1] * maxAmp, [0,1], 1);
        mrMin = linmap_(mrMin', [0, 1] * maxAmp, [0,1], 1);            
% end
[nEvt, nChans] = size(mrMin);
if isempty(maxPair), maxPair = nChans; end
[trX, trY] = deal(nan([nEvt, nChans, nChans], 'single'));
for chY = 1:nChans
    vrY1 = mrMin(:,chY);
    vlY1 = vrY1>0 & vrY1<1;
    for chX = 1:nChans
        if abs(chX-chY) > maxPair, continue; end
        if chY > chX
            vrX1 = mrMin(:,chX);
        else
            vrX1 = mrMax(:,chX);
        end
        viPlot1 = find(vrX1>0 & vrX1<1 & vlY1);
        trX(viPlot1,chY,chX) = vrX1(viPlot1) + chX - 1;
        trY(viPlot1,chY,chX) = vrY1(viPlot1) + chY - 1;
    end
end
% plot projection
viPlot = find(~isnan(trX) & ~isnan(trY));
vrX = trX(viPlot);  vrX=vrX(:);
vrY = trY(viPlot);  vrY=vrY(:);
tr_dim = size(trX);
end %func


%--------------------------------------------------------------------------
function plot_FigHist_(S0)

if nargin<1, S0 = get(0, 'UserData'); end
S_clu = S0.S_clu; P = S0.P;
[hFig, S_fig] = get_fig_cache_('FigHist');

nBins_hist = 50; % @TODO: put this in param file

vrX = logspace(0, 4, nBins_hist);
vrY1 = isi_hist_(S0.iCluCopy, vrX); 
vcTitle = sprintf('Cluster %d', S0.iCluCopy);

% draw
if isempty(S_fig) %first time the iCluPaste is always empty
    S_fig.hAx = axes_new_(hFig);
    S_fig.hPlot1 = stairs(S_fig.hAx, nan, nan, 'k'); 
    S_fig.hPlot2 = stairs(S_fig.hAx, nan, nan, 'r');     
    xlim(S_fig.hAx, [1 10000]); %in msec
    grid(S_fig.hAx, 'on');
    xlabel(S_fig.hAx, 'ISI (ms)');
    ylabel(S_fig.hAx, 'Prob. Density');
    set(S_fig.hAx, 'XScale', 'log');
end
update_plot_(S_fig.hPlot1, vrX, vrY1);
if ~isempty(S0.iCluPaste)
    vrY2 = isi_hist_(S0.iCluPaste, vrX);
    vcTitle = sprintf('Cluster %d (black) vs %d (red)', S0.iCluCopy, S0.iCluPaste);
    update_plot_(S_fig.hPlot2, vrX, vrY2);
else
    update_plot_(S_fig.hPlot2, nan, nan);
end
title_(S_fig.hAx, vcTitle);

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function vnHist = isi_hist_(iClu1, vrX)
P = get0_('P');
vrTime1 = double(clu_time_(iClu1)) / P.sRateHz;
vnHist = hist(diff(vrTime1)*1000, vrX);
vnHist(end)=0;
vnHist = vnHist ./ sum(vnHist);
end


%--------------------------------------------------------------------------
function plot_FigIsi_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
[hFig, S_fig] = get_fig_cache_('FigIsi');

[vrX1, vrY1] = get_returnMap_(S0.iCluCopy, P);                        
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    S_fig.hPlot1 = plot(S_fig.hAx, nan, nan, 'ko');
    S_fig.hPlot2 = plot(S_fig.hAx, nan, nan, 'ro');
    set(S_fig.hAx, 'XScale','log', 'YScale','log');   
    xlabel('ISI_{k} (ms)'); ylabel('ISI_{k+1} (ms)');
    axis(S_fig.hAx, [1 10000 1 10000]);
    grid(S_fig.hAx, 'on');
    % show refractory line
    line(get(S_fig.hAx,'XLim'), P.spkRefrac_ms*[1 1], 'Color', [1 0 0]);
    line(P.spkRefrac_ms*[1 1], get(S_fig.hAx,'YLim'), 'Color', [1 0 0]);
end  
update_plot_(S_fig.hPlot1, vrX1, vrY1);
if ~isempty(S0.iCluPaste)    
    [vrX2, vrY2] = get_returnMap_(S0.iCluPaste, P);
    update_plot_(S_fig.hPlot2, vrX2, vrY2);
else
    update_plot_(S_fig.hPlot2, nan, nan);
end

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [vrX, vrY] = get_returnMap_(iClu, P)
vrTime1 = double(clu_time_(iClu)) / P.sRateHz;
vrIsi1 = diff(vrTime1 * 1000); % in msec
vrX = vrIsi1(1:end-1);
vrY = vrIsi1(2:end);
viShow = randperm(numel(vrX), min(P.nShow, numel(vrX)));
vrX = vrX(viShow);
vrY = vrY(viShow);
end


%--------------------------------------------------------------------------
function [viTime_clu1, viSpk_clu1] = clu_time_(iClu1)
% returns time in sec
[P, S_clu, viTime_spk] = get0_('P', 'S_clu', 'viTime_spk');
viSpk_clu1 = S_clu.cviSpk_clu{iClu1};
viTime_clu1 = viTime_spk(S_clu.cviSpk_clu{iClu1});
end %func


%--------------------------------------------------------------------------
function plot_FigMap_(S0)
if nargin<1, S0 = get(0, 'UserData'); end 
P = S0.P; S_clu = S0.S_clu;
[hFig, S_fig] = get_fig_cache_('FigMap');

mrWav1 = S_clu.tmrWav_clu(:,:,S0.iCluCopy);
vrVpp = squeeze_(max(mrWav1) - min(mrWav1));
mrVpp = repmat(vrVpp(:)', [4, 1]);
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
    [S_fig.mrPatchX, S_fig.mrPatchY] = probe_map_(P);
    S_fig.hPatch = patch(S_fig.mrPatchX, S_fig.mrPatchY, mrVpp, ...
        'EdgeColor', 'k', 'FaceColor', 'flat');
    S_fig.alim = [min(S_fig.mrPatchX(:)), max(S_fig.mrPatchX(:)), min(S_fig.mrPatchY(:)), max(S_fig.mrPatchY(:))];
    colormap jet;
    mouse_figure(hFig);   
    nSites = size(P.mrSiteXY,1);    
    csText = arrayfun(@(i)sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
    S_fig.hText = text(P.mrSiteXY(:,1), P.mrSiteXY(:,2), csText, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');    
    xlabel('X Position (\mum)');
    ylabel('Y Position (\mum)');
else    
    set(S_fig.hPatch, 'CData', mrVpp);    
end
title_(S_fig.hAx, sprintf('Max: %0.1f \\muVpp', max(vrVpp)));
axis(S_fig.hAx, S_fig.alim);
caxis(S_fig.hAx, [0, max(vrVpp)]);

set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [mrPatchX, mrPatchY] = probe_map_(P)
vrX = [0 0 1 1] * P.vrSiteHW(2); 
vrY = [0 1 1 0] * P.vrSiteHW(1);
mrPatchX = bsxfun(@plus, P.mrSiteXY(:,1)', vrX(:));
mrPatchY = bsxfun(@plus, P.mrSiteXY(:,2)', vrY(:));
end %func


%--------------------------------------------------------------------------
function clu_info_(S0)
% This also plots cluster position
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
mh_info = findobj('Type', 'uimenu', 'Tag', 'mh_info');
S_clu1 = get_cluInfo_(S0.iCluCopy);
if ~isempty(S0.iCluPaste)
    S_clu2 = get_cluInfo_(S0.iCluPaste);
    vcLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', ...
        S0.iCluCopy, S_clu.csNote_clu{S0.iCluCopy}, ...
        S0.iCluPaste, S_clu.csNote_clu{S0.iCluPaste});
    set(mh_info, 'Label', vcLabel);
else
    S_clu2 = [];
    vcLabel = sprintf('Unit %d "%s"', S0.iCluCopy, S_clu.csNote_clu{S0.iCluCopy});
    set(mh_info, 'Label', vcLabel);
end
plot_FigPos_(S_clu1, S_clu2);
end %func


%--------------------------------------------------------------------------
function S_cluInfo = get_cluInfo_(iClu)
global trFet_spk %tnWav_spk tnWav_raw

% determine cluster position
if isempty(iClu), S_cluInfo=[]; return; end
[S0, P, S_clu] = get0_();

iSite1 = S_clu.viSite_clu(iClu);
viSite = P.miSites(:, iSite1);

% get cluster position
% [vrX1_spk, vrY1_spk, vrVpp_spk, mrVpp, trWav] = spikePos_(clu_time_(iClu), viSite, mrWav, P);
% viSpk1 = S_clu.cviSpk_clu{iClu};
% mrFet = squeeze_(trFet_spk(:,1,viSpk1));
% vrY1_spk = mrFet(end, viSpk1)'; 
% vrX1_spk = mrFet(end-1, viSpk1)'; 
% mrWav_clu = mean(trWav,3);
% xyPos = [median(vrX1_spk), median(vrY1_spk)];
xyPos = [S_clu.vrPosX_clu(iClu), S_clu.vrPosY_clu(iClu)];
vcPos = sprintf('Unit %d (x,y):(%0.1f, %0.1f)[pix]', iClu, xyPos/P.um_per_pix);
% if P.fWav_raw_show
%     mrWav_clu = S_clu.trWav_raw_clu(:,:,iClu);    
% else
%     mrWav_clu = S_clu.trWav_spk_clu(:,:,iClu);
% end
if P.fWav_raw_show
    mrWav_clu = S_clu.tmrWav_raw_clu(:,viSite,iClu);    
else
    mrWav_clu = S_clu.tmrWav_clu(:,viSite,iClu);    
end
trWav = trWav_clu_(iClu, P.nSpk_show * 1); 
if P.fWav_raw_show
    trWav = fft_lowpass_(trWav, get_set_(P, 'fc_spkwar_show', []), P.sRateHz);
end
% S_cluInfo = makeStruct_(xyPos, iClu, mrWav_clu, viSite, vrX1_spk, vrY1_spk, vcPos, trWav);
S_cluInfo = makeStruct_(xyPos, iClu, mrWav_clu, viSite, vcPos, trWav);
try
    S_cluInfo.l_ratio = S_clu.vrLRatio_clu(iClu);
    S_cluInfo.isi_ratio = S_clu.vrIsiRatio_clu(iClu);
    S_cluInfo.iso_dist = S_clu.vrIsoDist_clu(iClu);
    S_cluInfo.snr = S_clu.vrSnr_clu(iClu);
    S_cluInfo.uVmin = S_clu.vrVmin_uv_clu(iClu);
    S_cluInfo.uVpp = S_clu.vrVpp_uv_clu(iClu);
catch    
end
end %func


%--------------------------------------------------------------------------
function plot_FigPos_(S_clu1, S_clu2)
[hFig, S_fig] = get_fig_cache_('FigPos');
[S0, P, S_clu] = get0_();

% plot waveform in space
if isempty(S_fig)
    S_fig.hAx = axes_new_(hFig);
else
    cla(S_fig.hAx); hold(S_fig.hAx, 'on');
end
plot_unit_(S_clu1, S_fig.hAx, [0 0 0]);
%vrPosXY1 = [S_clu.vrPosX_clu(S_clu1.iClu), S_clu.vrPosY_clu(S_clu1.iClu)] / P.um_per_pix;
vrPosXY1 = [S_clu.vrPosX_clu(S_clu1.iClu), S_clu.vrPosY_clu(S_clu1.iClu)];
if isempty(S_clu2)    
    vcTitle = sprintf('Unit %d: (X=%0.1f, Y=%0.1f) [um]', S_clu1.iClu, vrPosXY1);
    try
        vcTitle = sprintf('%s\n%0.1fuVmin, %0.1fuVpp, SNR:%0.1f IsoDist:%0.1f ISIrat:%0.2f L-rat:%0.1f', ...
            vcTitle, S_clu1.uVmin, S_clu1.uVpp, S_clu1.snr, S_clu1.iso_dist, S_clu1.isi_ratio, S_clu1.l_ratio);
    catch
    end
else
    vrPosXY2 = [S_clu.vrPosX_clu(S_clu2.iClu), S_clu.vrPosY_clu(S_clu2.iClu)] / P.um_per_pix;
    plot_unit_(S_clu2, S_fig.hAx, [1 0 0]);
    vcTitle = sprintf('Unit %d/%d(red)\n(X=%0.1f/%0.1f, Y=%0.1f/%0.1f) [um]', ...
        S_clu1.iClu, S_clu2.iClu, ...
        [vrPosXY1(1), vrPosXY2(1), vrPosXY1(2), vrPosXY2(2)]);
end
title_(S_fig.hAx, vcTitle);
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function plot_unit_(S_clu1, hAx, vcColor0)
if isempty(S_clu1), return; end
if nargin<2, hAx = axes_new_('FigWav'); end
if nargin<3, vcColor0 = [0 0 0]; end
[S0, P, S_clu] = get0_();
[~, S_figWav] = get_fig_cache_('FigWav');
maxAmp = S_figWav.maxAmp;
% plot individual unit
nSamples = size(S_clu1.mrWav_clu,1);
vrX = (1:nSamples)'/nSamples;
vrX([1,end])=nan; % line break

if ~isequal(vcColor0, [0 0 0])
    trWav1 = zeros(1,1,0);
else
    trWav1 = S_clu1.trWav;
end

% show example traces
for iWav = size(trWav1,3):-1:0
    if iWav==0
        mrY1 = S_clu1.mrWav_clu / maxAmp;
        lineWidth=1.5;
        vcColor = vcColor0;
    else
        mrY1 = trWav1(:,:,iWav) / maxAmp;
        lineWidth=.5;
        vcColor = .5*[1,1,1];
    end
    vrX1_site = P.mrSiteXY(S_clu1.viSite, 1) / P.um_per_pix;
    vrY1_site = P.mrSiteXY(S_clu1.viSite, 2) / P.um_per_pix;
    mrY1 = bsxfun(@plus, mrY1, vrY1_site');
    mrX1 = bsxfun(@plus, repmat(vrX, [1, size(mrY1, 2)]), vrX1_site');
    line(mrX1(:), mrY1(:), 'Color', vcColor, 'Parent', hAx, 'LineWidth', lineWidth);
end
xlabel(hAx, 'X pos [pix]');
ylabel(hAx, 'Z pos [pix]');
grid(hAx, 'on');
xlim(hAx, [min(mrX1(:)), max(mrX1(:))]);
ylim(hAx, [floor(min(mrY1(:))-1), ceil(max(mrY1(:))+1)]);
end %func


%--------------------------------------------------------------------------
function trWav1 = trWav_clu_(iClu1, nSpk_show)
% Get a subset of spike waveforms of a cluster
global tnWav_spk tnWav_raw

if nargin<2, nSpk_show=inf; end

S0 = get(0, 'UserData');
P = S0.P;
[viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(S0.S_clu, iClu1, S0.viSite_spk);
viSpk_clu1 = randomSelect_(viSpk_clu1, nSpk_show);
if P.fWav_raw_show
    trWav1 = raw2uV_(tnWav_raw(:,:,viSpk_clu1), P);
else
    trWav1 = tnWav2uV_(tnWav_spk(:,:,viSpk_clu1), P);
end
end %func


%--------------------------------------------------------------------------
function rescale_fig_(event, S0, P)
figure_wait_(1);
set(0, 'UserData', S0);

switch (get(gcf, 'Tag'))
    case 'FigWav'  % rescale all the figures
        [S_fig, maxAmp_prev, hFigWav] = set_fig_maxAmp_('FigWav', event);                
%         delete(S_fig.vhPlot);
%         rescale_FigWav_(S0.S_clu, P);
        set_fig_(hFigWav, plot_tnWav_clu_(S_fig, P));
        multiplot(S0.hCopy, S_fig.maxAmp);
        if ~isempty(S0.iCluPaste)
            multiplot(S0.hPaste, S_fig.maxAmp);
        end
        rescale_spikes_(S_fig.hSpkAll, maxAmp_prev, P);
        title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp)); %update scale

    case 'FigTime' % set time view        
        [S_fig, maxAmp_prev] = set_fig_maxAmp_('FigTime', event);
        switch lower(P.vcFet_show)    
            case 'vpp'
                ylim(S_fig.hAx, [0 S_fig.maxAmp]);
                imrect_set_(S_fig.hRect, [], [0 S_fig.maxAmp*2]);
            otherwise
                ylim(S_fig.hAx, [0, 1] * S_fig.maxAmp);
                imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);
        end

    case 'FigProj' % set projection view        
        [S_fig, maxAmp_prev] = set_fig_maxAmp_('FigProj', event);
        keyPressFcn_FigProj_([], event);
end %switch
figure_wait_(0);
end


%--------------------------------------------------------------------------
function button_CluWav_(xyPos, vcButton)
if strcmpi(vcButton, 'normal')
    event.Button = 1;
elseif strcmpi(vcButton, 'alt')
    event.Button = 3;
else
    return;
end
xPos = round(xyPos(1));
S0 = get(0, 'UserData');
switch(event.Button)
    case 1 %left click. copy clu and delete existing one
        S0 = update_cursor_(S0, xPos, 0);
    case 2 %middle, ignore
        return; 
    case 3 %right click. paste clu
        S0 = update_cursor_(S0, xPos, 1);
end
figure_wait_(1);
set0_(S0);
keyPressFcn_cell_(get_fig_cache_('FigWav'), {'c','t','j','i','v','e','f'}, S0); %'z'
plot_raster_();
auto_scale_proj_time_();
figure_wait_(0);
end %func


%--------------------------------------------------------------------------
% function rescale_FigWav_(S_clu, P)
% S0 = get(0, 'UserData');
% [hFig, S_fig] = get_fig_cache_('FigWav');
% vhPlot = S_fig.vhPlot;
% % tmrCluWav = hideCluSite_(Sclu, P);
% for iPlot=1:numel(vhPlot)
%     
% end %for
% % for iClu=1:S_clu.nClu
% %     viSites1 = P.miSites(:, S_clu.viSite_clu(iClu));
% %     mrY1 = S_clu.tmrWav_clu(:,viSites1,iClu) / S_fig.maxAmp; 
% %     mrY1 = bsxfun(@plus, mrY1, single(viSites1'));
% %     set(vhPlot(iClu), 'YData', mrY1(:));
% % end
% end %func


%--------------------------------------------------------------------------
function keyPressFcn_FigTime_(hObject, event, S0)
% global mrWav
if nargin<3, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
hFig = hObject;
S_fig = get(hFig, 'UserData');
% mrWav = S.mrWav;
nSites = numel(P.viSite2Chan);
% set(hObject, 'Pointer', 'watch');
% figure_wait(1);
switch lower(event.Key)
    case {'leftarrow', 'rightarrow'}
        if ~isVisible_(S_fig.hAx)
            msgbox_('Channel switching is disabled in the position view'); return; 
        end
        factor = key_modifier_(event, 'shift')*3 + 1;
        if strcmpi(event.Key, 'rightarrow')
            S_fig.iSite = min(S_fig.iSite + factor, nSites);
        else
            S_fig.iSite = max(S_fig.iSite - factor, 1);
        end
        set(hFig, 'UserData', S_fig);        
        update_FigTime_();                

    case {'uparrow', 'downarrow'} %change ampl
        if ~isVisible_(S_fig.hAx)
            msgbox_('Zoom is disabled in the position view'); return; 
        end
        rescale_fig_(event, S0, P);
        
    case 'r' %reset view
        if ~isVisible_(S_fig.hAx), return; end
        axis(S_fig.hAx, [S_fig.time_lim, S_fig.vpp_lim]);
        imrect_set_(S_fig.hRect, S_fig.time_lim, S_fig.vpp_lim);

    case 'm' %merge
        ui_merge_();
        
    case 'h' %help
        msgbox_(S_fig.csHelp, 1);
        
    case 'b' %background spike toggle
        if isVisible_(S_fig.hAx)  
            S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0);
        else
            S_fig.fPlot0 = toggleVisible_(S_fig.hPlot0_track);
        end
        set(hFig, 'UserData', S_fig);
        
    case 't'
        plot_FigTime_(S0);
        
    case 'z' % track depth
        disp('FigTime:''z'' not implemented yet');
%         plot_SpikePos_(S0, event);
        
    case 's' %split. draw a polygon
        if ~isempty(S0.iCluPaste)
            msgbox_('Select one cluster'); return;
        end
        try
            hPoly = impoly_();
            if isempty(hPoly); return ;end
            mrPolyPos = getPosition(hPoly);
            vrX1 = double(get(S_fig.hPlot1, 'XData'));
            vrY1 = double(get(S_fig.hPlot1, 'YData'));
            vlIn = inpolygon(vrX1, vrY1, mrPolyPos(:,1), mrPolyPos(:,2));
            hSplit = line(vrX1(vlIn), vrY1(vlIn), 'Color', [1 0 0], 'Marker', '.', 'LineStyle', 'none');
            if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'yes')
                split_clu_(S0.iCluCopy, vlIn);
            end
            delete_multi_(hPoly, hSplit);
        catch
            disp(lasterror());
        end
        
    case 'p' %update projection view
        vrPos = getPosition(S_fig.hRect);
        tlim_proj = [vrPos(1), sum(vrPos([1,3]))];
        P.tlim_proj = tlim_proj;
        plot_FigProj_(S0);
        
%     case 'f' % feature display instead of amplitude display
%         if strcmpi(P.vcFet_show, 'fet')
%             P.vcFet_show = 'vpp';
%         else
%             P.vcFet_show = 'fet';
%         end
%         set0_(P);
%         update_FigTime_();     
        
    case 'c' % compare pca across channels
        disp('FigTime: Not implemented yet'); return;
%         hMsg = msgbox_('Plotting...');
%         figure; hold on;
%         [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(S0.iCluCopy);
%         [~, mrPv1] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);
%         mrPv1 = norm_mr_(mrPv1);
%         
%         if key_modifier_(event, 'control') %show chain of clusters
%             trPv1 = mrPv1;
%             iClu_next = get_next_clu_(S_clu, S0.iCluCopy);
%             viClu_track = S0.iCluCopy;
%             while ~isempty(iClu_next)
%                 [mrWav_mean1, viSite1] = mrWav_int_mean_clu_(iClu_next);
%                 [~, mrPv1a] = pca(mrWav_mean1, 'NumComponents', P.nPc_dip, 'Center', 1);        
%                 mrPv1a = norm_mr_(mrPv1a);
%                 mrPv1 = flip_prinvec_(mrPv1a, mean(trPv1,3));
%                 trPv1 = cat(3, trPv1, mrPv1);
%                 viClu_track(end+1) = iClu_next;
%                 
%                 iClu_next = get_next_clu_(S_clu, iClu_next);
%             end      
%             multiplot(plot(nan,nan,'k'), 1, 1:size(trPv1,1), trPv1);
% %             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
%             vcTitle = sprintf('PCA across chan: Clu %s', sprintf('%d,', viClu_track));
%         elseif ~isempty(S0.iCluPaste)            
%             [mrWav_mean2, viSite1] = mrWav_int_mean_clu_(S0.iCluPaste);
%             [~, mrPv2] = pca(mrWav_mean2, 'NumComponents', P.nPc_dip);     
%             mrPv2 = match_mrPv_(mrPv2, mrPv1);
% %             mrPv2 = flip_prinvec_(mrPv2, mrPv1);
%             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'k');
%             mr2plot(norm_mr_(mrPv2), 'scale', 1, 'LineStyle', 'r--');            
%             vcTitle = sprintf('PCA across chan: Clu %d vs %d', S0.iCluCopy, S0.iCluPaste);            
%         else        
%             mr2plot(norm_mr_(mrPv1), 'scale', 1, 'LineStyle', 'r');
%             vcTitle = sprintf('PCA across chan: Clu %d', S0.iCluCopy);                        
%         end
% %         mr2plot(mrPv1, 'scale', 1, 'LineStyle', 'k');
%         grid on; 
%         title_(vcTitle);
% %         if ~isempty(S0.iCluPaste)   
% %             compare_interp_(Sclu, S0.iCluCopy, S0.iCluPaste);
% %         end
%         try close(hMsg); catch; end
        
%     case 'f' %feature export
%         eval(sprintf('mrFet_clu%d = getFet_clu_(S0.iCluCopy);', S0.iCluCopy));
%         mrDist1 = squareform(pdist(mrFet1'));
%         vrFet1 = sqrt(sum(mrFet1.^2));
%         mrDist1 = bsxfun(@rdivide, mrDist1, vrFet1); %norm        
%         eval(sprintf('assignWorkspace_(mrFet_clu%d);', S0.iCluCopy));
        
    case 'e' %export selected to workspace
        disp('FigTime: ''e'' not implemented yet'); return;
end %switch
% drawnow;
% figure_wait(0);
end %func


%--------------------------------------------------------------------------
function rescale_spikes_(hSpkAll, maxAmp_prev, P)
S = get(hSpkAll, 'UserData');
[~, S_fig] = get_fig_cache_('FigWav');
S0 = get(0, 'UserData');
cvrY = S.cvrY;
cviSite = S.cviSite;
% nSamples = diff(P.spkLim)+1;
scale = S_fig.maxAmp / maxAmp_prev;
for iClu=1:numel(cvrY)
    viSite1 = cviSite{iClu};
    nSites1 = numel(viSite1);
    trY = reshape(cvrY{iClu}, [], nSites1, S.vnSpk(iClu));
    for iSite1 = 1:nSites1
        y_off = viSite1(iSite1);
        trY(:,iSite1,:) = (trY(:,iSite1,:) - y_off) / scale + y_off;
    end
    cvrY{iClu} = trY(:);
end
S.cvrY = cvrY;
set(hSpkAll, 'YData', cell2mat_(cvrY), 'UserData', S);
end


%--------------------------------------------------------------------------
function keyPressFcn_FigProj_(hFig, event)
% hFig = hObject;
% global mrWav
S0 = get(0, 'UserData');
[P, S_clu] = get0_('P', 'S_clu');
[hFig, S_fig] = get_fig_cache_('FigProj');
% nSites = numel(P.viSite2Chan);
S_plot1 = get(S_fig.hPlot1, 'UserData');
viSites_show = S_plot1.viSites_show;
% nSites = numel(viSites_show);
% set(hObject, 'Pointer', 'watch');
figure_wait_(1);
switch lower(event.Key)
    case {'uparrow', 'downarrow'}
        rescale_FigProj_(hFig, event, S_fig);

    case {'leftarrow', 'rightarrow'} % change channels
        fPlot = 0;
        if strcmpi(event.Key, 'leftarrow')
            if min(S_fig.viSites_show)>1
                S_fig.viSites_show=S_fig.viSites_show-1; 
                fPlot = 1;
            end
        else
            if max(S_fig.viSites_show) < max(P.viSite2Chan)
                S_fig.viSites_show=S_fig.viSites_show+1;                 
                fPlot = 1;
            end
        end
        if fPlot
            set(hFig, 'UserData', S_fig);
            S0.P.viSites_show = S_fig.viSites_show;
            plot_FigProj_(S0);
        end
        
    case 'r' %reset view
        axis([0 numel(viSites_show) 0 numel(viSites_show)]);

    case 's' %split
        figure_wait_(0);
        if ~isempty(S0.iCluPaste)
            msgbox_('Select one cluster to split'); return;
        end
        S_plot1 = select_polygon_(S_fig.hPlot1); 
        if ~isempty(S_plot1)
            [fSplit, vlIn] = plot_split_(S_plot1);
            if fSplit
                S_clu = split_clu_(S0.iCluCopy, vlIn);
            else
                update_plot2_proj_();
%                 delete_multi_(S_plot1.hPoly);
            end
        end
        
    case 'm'
        ui_merge_();
        
    case 'f'
        disp('keyPressFcn_FigProj_: ''f'': not implemented yet');
%         if strcmpi(P.vcFet_show, 'vpp')
%             S0.vcFet_show = P.vcFet;
%         else
%             S0.vcFet_show = 'vpp';
%         end
%         set(0, 'UserData', S0);        
%         plot_FigProj_();        
        
    case 'b' %background spikes
        toggleVisible_(S_fig.hPlot0);

    case 'h' %help
        msgbox_(S_fig.csHelp, 1);
end %switch
% drawnow;
figure_wait_(0);
end %func


%--------------------------------------------------------------------------
function S_fig = rescale_FigProj_(hFig, event, S_fig)
% S_fig = rescale_FigProj_(hFig, event)
% S_fig = rescale_FigProj_(hFig, event, S_fig)

S0 = get(0, 'UserData');
if nargin<3, S_fig = get(hFig, 'UserData'); end
S_fig.maxAmp = change_amp_(event, S_fig.maxAmp);     
vhPlot = [S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2];
if isempty(S0.iCluPaste), vhPlot(end) = []; end
rescaleProj_(vhPlot, S_fig.maxAmp);
xlabel(sprintf(S_fig.vcXLabel, S_fig.maxAmp));   
ylabel(sprintf(S_fig.vcYLabel, S_fig.maxAmp));  
if nargout==0
    set(hFig, 'UserData', S_fig);
end
end


%--------------------------------------------------------------------------
function rescaleProj_(vhPlot1, maxAmp)
P = get0_('P');
for iPlot=1:numel(vhPlot1)
    hPlot1 = vhPlot1(iPlot);
    S_plot1 = get(hPlot1, 'UserData');
    update_plot2_proj_();
    S_plot1 = struct_delete_(S_plot1, 'hPoly'); %, 'hPlot_split'
    [vrX, vrY, viPlot, tr_dim] = amp2proj_(S_plot1.mrMin, S_plot1.mrMax, maxAmp, P.maxSite_show);
    S_plot1 = struct_add_(S_plot1, viPlot, vrX, vrY, maxAmp);
    set(hPlot1, 'XData', vrX, 'YData', vrY, 'UserData', S_plot1);
end
end %func


%--------------------------------------------------------------------------
function button_FigWavCor_(xyPos, vcButton)
S0 = get(0, 'UserData');
xyPos = round(xyPos);
switch lower(vcButton)
    case 'normal' %left click        
        S0.iCluCopy = xyPos(1);
        if diff(xyPos) == 0
            S0.iCluPaste = [];
        else
            S0.iCluPaste = xyPos(2);
        end
        S0 = button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
        S0 = keyPressFcn_cell_(get_fig_cache_('FigWav'), {'z'}, S0); %zoom        
end %switch
end %func


%--------------------------------------------------------------------------
function [mrFet1, viSpk1] = getFet_clu_(iClu1, iSite)
% get features on-fly
% global tnWav_spk
% S0 = get(0, 'UserData');
MAX_SAMPLE = 40000;     % max points ti display
[S_clu, P, viSite_spk] = get0_('S_clu', 'P', 'viSite_spk');
% if nargin<2, viSite = P.miSites(:, S0.S_clu.viSite_clu(iClu1)); end
if isempty(iClu1) % select spikes based on sites
    n_use = 1 + round(P.maxSite);
    viSpk1 = find(ismember(viSite_spk, P.miSites(1:n_use, iSite)));
    viSpk1 = randomSelect_(viSpk1, MAX_SAMPLE);    
else
    viSpk1 = S_clu.cviSpk_clu{iClu1};
end

switch lower(P.vcFet_show)
    case {'vmin', 'vpp'}
        mrWav_spk1 = squeeze_(tnWav2uV_(tnWav_spk_sites_(viSpk1, iSite)));
        mrFet1 = max(mrWav_spk1)-min(mrWav_spk1);
    case 'cov'
        mrFet1 = calc_cov_spk_(viSpk1, iSite);
    case 'pca'
        mrFet1 = pca_pc_spk_(viSpk1, iSite);
    otherwise
        error('not implemented yet');
end
mrFet1 = squeeze_(abs(mrFet1));
end %func


%--------------------------------------------------------------------------
function [vi1, vrX1, vrY1, xlim1, ylim1] = imrect_plot_(hRect, hPlot1)
% return points inside
vrX = get(hPlot1, 'XData');
vrY = get(hPlot1, 'YData');
vrPos_rect = getPosition(hRect);
xlim1 = vrPos_rect(1) + [0, vrPos_rect(3)];
ylim1 = vrPos_rect(2) + [0, vrPos_rect(4)];
vi1 = find(vrX>=xlim1(1) & vrX<=xlim1(2) & vrY>=ylim1(1) & vrY<=ylim1(2));
if nargout>=2
    vrX1 = vrX(vi1);
    vrY1 = vrY(vi1);
end
end %func


%--------------------------------------------------------------------------
function [tnWav_raw, tnWav_spk] = mn2tn_wav_(mnWav_raw, mnWav_spk, viSite_spk, viTime_spk, P)
nSpks = numel(viTime_spk);
nSites = numel(P.viSite2Chan);
spkLim_raw = P.spkLim * 2;
spkLim_wav = P.spkLim;
nSites_spk = (P.maxSite * 2) + 1;
tnWav_raw = zeros(diff(spkLim_raw) + 1, nSites_spk, nSpks, 'int16');
tnWav_spk = zeros(diff(spkLim_wav) + 1, nSites_spk, nSpks, 'int16');
if isempty(viSite_spk)
    tnWav_raw = permute(gather_(mr2tr3_(mnWav_raw, spkLim_raw, viTime_spk)), [1,3,2]);
    tnWav_spk = permute(gather_(mr2tr3_(mnWav_spk, spkLim_wav, viTime_spk)), [1,3,2]);
else
    for iSite = 1:nSites
        viiSpk11 = find(viSite_spk == iSite);
        if isempty(viiSpk11), continue; end
        viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
        viSite11 = P.miSites(:,iSite);
        try
            tnWav_raw(:,:,viiSpk11) = permute(gather_(mr2tr3_(mnWav_raw, spkLim_raw, viTime_spk11, viSite11)), [1,3,2]); %raw
            tnWav_spk(:,:,viiSpk11) = permute(gather_(mr2tr3_(mnWav_spk, spkLim_wav, viTime_spk11, viSite11)), [1,3,2]);
        catch
            mnWav_raw = gather_(mnWav_raw);
            mnWav_spk = gather_(mnWav_spk);
            tnWav_raw(:,:,viiSpk11) = permute(mr2tr3_(mnWav_raw, spkLim_raw, viTime_spk11, viSite11), [1,3,2]); %raw
            tnWav_spk(:,:,viiSpk11) = permute(mr2tr3_(mnWav_spk, spkLim_wav, viTime_spk11, viSite11), [1,3,2]);
        end
    end
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Bugfix: returning S0 so that the cluster can be updated
function S0 = ui_delete_()
S0 = get(0, 'UserData');
if ~isempty(S0.iCluPaste)
    msgbox_('Must select one cluster'); return;
end        
figure_wait_(1);

iClu_del = S0.iCluCopy;
hMsg = msgbox_open_('Deleting...');
S0.S_clu = delete_clu_(S0.S_clu, S0.iCluCopy);
set(0, 'UserData', S0);
plot_FigWav_(S0); %redraw plot 
% S0.S_clu.mrWavCor = wavCor_delete_(S0.iCluCopy); 
FigWavCor_update_(S0);
S0.iCluCopy = min(S0.iCluCopy, S0.S_clu.nClu);
set(0, 'UserData', S0);
button_CluWav_simulate_(S0.iCluCopy);

close_(hMsg);
figure_wait_(0);
fprintf('%s [W] deleted Clu %d\n', datestr(now, 'HH:MM:SS'), iClu_del);
S0 = save_log_(sprintf('delete %d', iClu_del), S0);
set(0, 'UserData', S0);
end


%--------------------------------------------------------------------------
function FigWavCor_update_(S0)
if nargin<1, S0 = get0_(); end

[hFig, S_fig] = get_fig_cache_('FigWavCor');
set(S_fig.hImWavCor, 'CData', S0.S_clu.mrWavCor);
% plotDiag_([0, nSites], '-', 'Color', [0 0 0], 'LineWidth', 1.5, 'Parent', S_fig.); %plot in one scoop
nClu = size(S0.S_clu.mrWavCor, 1);
[vrX, vrY] = plotDiag__([0, nClu, .5]);
set(S_fig.hDiag, 'XData', vrX, 'YData', vrY);
end %func


%--------------------------------------------------------------------------
function S_clu = delete_clu_(S_clu, iClu)
% sets the cluster to zero
% iClu0 = Sclu.viCluNum(iClu);
% S_clu = struct_set_(S_clu, iClu, [], ...
%     'cviSpk_clu', 'vnSpk_clu', 'viSite_clu', 'csNote_clu', 'vrPosX_clu', 'vrPosY_clu');

vlKeep_clu = true(S_clu.nClu, 1);
vlKeep_clu(iClu) = 0;
S_clu = S_clu_select_(S_clu, vlKeep_clu);
S_clu.mrWavCor = S_clu.mrWavCor(vlKeep_clu, vlKeep_clu);
% S_clu.tmrWav_clu(:,:,iClu) = [];
% try S_clu.tmrWav_raw_clu(:,:,iClu) = []; catch; end

iClu_del = min(S_clu.viClu) - 1;
if iClu_del==0, iClu_del = -1; end
S_clu.viClu(S_clu.viClu==iClu) = iClu_del;
S_clu.nClu = numel(S_clu.vnSpk_clu);
% update viClu
if iClu < max(S_clu.viClu)
    viUpdate = find(S_clu.viClu>iClu);
    S_clu.viClu(viUpdate) = S_clu.viClu(viUpdate) - 1;
end
for iClu3 = iClu+1:S_clu.nClu % update cluster chain info
    S_clu = S_clu_update_note_(S_clu, iClu3, get_next_clu_(S_clu, iClu3) - 1);
end
assert(S_clu_valid_(S_clu), 'Cluster number is inconsistent after deleting');
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_update_note_(S_clu, iClu1, iClu_next)
if isempty(iClu_next), return ;end
vcNote_clu1 = S_clu.csNote_clu{iClu1};
if isempty(vcNote_clu1), return; end
iStart = find(vcNote_clu1 == '=', 1, 'first');
if isempty(iStart), return; end
vcPre = vcNote_clu1(1:iStart);
vcNote_clu1 = vcNote_clu1(iStart+1:end);
iEnd = find(vcNote_clu1 == ',' | vcNote_clu1 == ';' | vcNote_clu1 == ' ', 1, 'first');
if ~isempty(iEnd)
    vcNote_clu1 = vcNote_clu1(1:iEnd-1); 
    vcPost = vcNote_clu1(iEnd:end);
else
    vcPost = '';
end
S_clu.csNote_clu{iClu1} = sprintf('%s%d%s', vcPre, iClu_next, vcPost);

if isnan(iClu_next), iClu_next = []; return; end
end %func


%--------------------------------------------------------------------------
function iClu_next = get_next_clu_(S_clu, iClu1)
%  get a next cluster number reading from Sclu.csNote_clu "=Clu#"
iClu_next = [];
vcNote_clu1 = S_clu.csNote_clu{iClu1};
if isempty(vcNote_clu1), return; end
iStart = find(vcNote_clu1 == '=', 1, 'first');
if isempty(iStart), return; end
vcNote_clu1 = vcNote_clu1(iStart+1:end);
iEnd = find(vcNote_clu1 == ',' | vcNote_clu1 == ';' | vcNote_clu1 == ' ', 1, 'first');
if ~isempty(iEnd), vcNote_clu1 = vcNote_clu1(1:iEnd-1); end
iClu_next = str2double(vcNote_clu1);
if isnan(iClu_next), iClu_next = []; return; end
end %func


%--------------------------------------------------------------------------
function restore_clu_(varargin)
% global mrWav
% restore last deleted. most negative clu is last deleted
% error('to be fixed. Fix centroid code');
[S0, P, S_clu] = get0_();
iClu_del = min(S_clu.viClu);
if iClu_del >=0, msgbox_('Deleted cluster is not found'); return; end

figure_wait_(1);
% if deleted add a clu at the end and zoom at it
% change clusters
iClu_new = double(max(S_clu.viClu) + 1);
S_clu.viClu(S_clu.viClu == iClu_del) = iClu_new;
S_clu.nClu = iClu_new;
S_clu = S_clu_update_(S_clu, iClu_new, P);
% S_clu = S_clu_position_(S_clu, iClu_new);
S_clu.csNote_clu{end+1} = '';
[S_clu, iClu_new] = clu_reorder_(S_clu);

% update all the other views
% delete_multi_(S0.vhPlot, S0.vhText);
% S0.S_clu = S_clu; set(0, 'UserData', S0);
[S_clu, S0] = S_clu_commit_(S_clu);
plot_FigWav_(S0); %redraw plot
plot_FigWavCor_(S0);
% S0 = set0_(mrWavCor);
set(0, 'UserData', S0);

% append to the end for now
button_CluWav_simulate_(iClu_new);
keyPressFcn_cell_(get_fig_cache_('FigWav'), 'z');
fprintf('%s [W] Restored Clu %d\n', datestr(now, 'HH:MM:SS'), iClu_new);
figure_wait_(0);
end


%--------------------------------------------------------------------------
function [vrY, vrX] = tr2plot_(trWav, iPos, viSite_show, maxAmp, P)
if nargin<2, iPos=1; end
iPos = double(iPos);
if nargin<5, P = get0_('P'); end %S0 = get(0, 'UserData'); P = S0.P;
% [~, S_fig] = get_fig_cache_('FigWav');
% if isfield(S_fig, 'maxAmp')
%     maxAmp = S_fig.maxAmp;    
% else
%     maxAmp = P.maxAmp;
% end

if nargin<3, viSite_show = []; end
% P = funcDefStr_(P, 'LineStyle', 'k', 'spkLim', [-10 24], 'maxAmp', 500, 'viSite_show', []);
% P.LineStyle
% if isempty(P.LineStyle), P.LineStyle='k'; end
if isempty(viSite_show), viSite_show = 1:size(trWav,2); end

[nSamples, nChans, nSpk] = size(trWav);
nSites_show = numel(viSite_show);
trWav = single(trWav) / maxAmp;
trWav = trWav + repmat(single(viSite_show(:)'), [size(trWav,1),1,size(trWav,3)]);
trWav([1,end],:,:) = nan;
vrY = trWav(:);

if nargout>=2
    x_offset = (P.spkLim(2))/(diff(P.spkLim)+1) + iPos - 1;
    vrX = (1:nSamples)/nSamples + x_offset;
    vrX([1,end]) = nan;
    vrX = repmat(vrX, [1, nSites_show * nSpk]);
    vrX = single(vrX(:));
end
end


%--------------------------------------------------------------------------
% 8/7/17 JJJ: tested and documented, file creation date
function [csFile_merge, vcDir] = dir_file_(vcFile_dir, fSortByDate)
% function [csFile_merge, vcDir, csFile_merge1] = dir_file_(vcFile_dir, fSortByDate)
% search for files and sort by date

if nargin<2, fSortByDate=1; end

[vcDir, ~, ~] = fileparts(vcFile_dir);
if isempty(vcDir), vcDir = '.'; end
csFile_merge = dir_(vcFile_dir);

% vsDir = dir(vcFile_dir);
% vrDatenum = cell2mat_({vsDir.datenum});
% csFile_merge = {vsDir.name};
switch fSortByDate
    case 0
        ; %no change
    case 1
        vrDatenum = file_created_meta_(csFile_merge, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);        
    case 2
        vrDatenum = file_created_(csFile_merge, 'creationTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);
    case 3
        vrDatenum = file_created_(csFile_merge, 'lastModifiedTime');
        [~,ix] = sort(vrDatenum, 'ascend');
        csFile_merge = csFile_merge(ix);
    case 4
        [csFile_merge, ix] = sort_nat_(csFile_merge, 'ascend');        
    otherwise
        fprintf(2, 'dir_file_: Invalid option: %d\n', fSortByDate);
end
% csFile_merge1 = csFile_merge; 
% csFile_merge = cellfun(@(vc)[vcDir, filesep(), vc], csFile_merge, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function S0 = ui_merge_()
[S0, P] = get0_();

if isempty(S0.iCluPaste)
    msgbox_('Right-click a cluster to merge.', 1); return;
end
if S0.iCluCopy == S0.iCluPaste
    msgbox_('Cannot merge to itself.', 1); return;
end

figure_wait_(1);
S0.S_clu = merge_clu_(S0.S_clu, S0.iCluCopy, S0.iCluPaste, P);
set(0, 'UserData', S0);
plot_FigWav_(S0); %redraw plot
S0.iCluCopy = min(S0.iCluCopy, S0.iCluPaste);
S0.iCluPaste = [];
set(0, 'UserData', S0);
update_plot_(S0.hPaste, nan, nan);
S0 = update_FigCor_(S0);        
S0 = button_CluWav_simulate_(S0.iCluCopy, [], S0);
S0 = save_log_(sprintf('merge %d %d', S0.iCluCopy, S0.iCluPaste), S0);
set(0, 'UserData', S0);

% msgbox_close(hMsg);
figure_wait_(0);
% S_clu = S0.S_clu;
end %func


%--------------------------------------------------------------------------
function S_clu = merge_clu_(S_clu, iClu1, iClu2, P)
if iClu1>iClu2, [iClu1, iClu2] = swap_(iClu1, iClu2); end

S_clu = merge_clu_pair_(S_clu, iClu1, iClu2);
S_clu = S_clu_refrac_(S_clu, P, iClu1); % remove refrac
S_clu = S_clu_update_(S_clu, iClu1, P);
S_clu = delete_clu_(S_clu, iClu2);
% S_clu = S_clu_remove_empty_(S_clu);
assert(S_clu_valid_(S_clu), 'Cluster number is inconsistent after merging');
fprintf('%s [W] merging Clu %d and %d\n', datestr(now, 'HH:MM:SS'), iClu1, iClu2);
end %func


%--------------------------------------------------------------------------
function S0 = update_FigCor_(S0)
if nargin<1, S0 = get(0, 'UserData'); end
P = S0.P; S_clu = S0.S_clu;
[hFig, S_fig] = get_fig_cache_('FigWavCor');

% figure(hFig);   
% hMsg = msgbox_open('Computing Correlation');

xylim = get(gca, {'XLim', 'YLim'});
S0.S_clu = S_clu;
plot_FigWavCor_(S0);
set(gca, {'XLim', 'YLim'}, xylim);

set(S_fig.hCursorV, 'XData', S0.iCluCopy*[1 1]);
set(S_fig.hCursorH, 'YData', S0.iCluCopy*[1 1]);

if nargout==0
    set(0, 'UserData', S0); %update field
end
end %func


%--------------------------------------------------------------------------
function vlSuccess = download_files_(csLink, csDest)
% download file from the web
if nargin<2, csDest = link2file_(csLink); end
vlSuccess = false(size(csLink));
for iFile=1:numel(csLink)    
    try
        % download from list of files    
        fprintf('\tDownloading %s: ', csLink{iFile});
        vcFile_out1 = websave(csDest{iFile}, csLink{iFile});
        fprintf('saved to %s\n', vcFile_out1);
        vlSuccess(iFile) = 1;
    catch
        fprintf(2, '\n\tCannot download. Check internet connection.\n');
    end
end %for
end %func


%--------------------------------------------------------------------------
% function auto_split_(fMulti)
% % Auto-split feature that calls Hidehiko Inagaki's code
% % 20160426
% % global tnWav_spk
% % error('auto_split_ @TODO mrWav dependency');
% if nargin<1, fMulti = 0; end
% 
% [S0, P, S_clu] = get0_(); %S0 = get(0, 'UserData');
% iClu1 = S0.iCluCopy;
% iSite1 = S_clu.viSite_clu(iClu1);
% 
% 
% % mrSpkWav1 = vr2mr2(mrWav, viTime1, S0.P.spkLim, iSite1);
% % mrSpkWav1 = tnWav2uV_(squeeze_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, iSite1)));
% mrSpkWav1 = tnWav2uV_(squeeze_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, iSite1)));
% 
% % preview clusters to split
% nSplit = preview_split_(mrSpkWav1);
% if isnan(nSplit), return; end
% vlSpkIn = auto_split_wav_(mrSpkWav1, nSplit);
% 
% hFigTemp = gcf;
% if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'yes')
%     split_clu_(iClu1, vlSpkIn);    
% end
% close(hFigTemp);
% end %func


%--------------------------------------------------------------------------
function auto_split_(fMulti)
% Auto-split feature that calls Hidehiko Inagaki's code
% 20160426
% global tnWav_spk
if nargin<1, fMulti = 0; end

S0 = get(0, 'UserData');
S_clu = S0.S_clu;
if ~isempty(S0.iCluPaste), msgbox_('Select one cluster'); return; end
if S_clu.vnSpk_clu(S0.iCluCopy)<2, msgbox_('At least two spikes required for splitting'); return; end
    
hMsg = msgbox_('splitting...');
iClu1 = S0.iCluCopy;
iSite1 = S0.S_clu.viSite_clu(iClu1);
if fMulti
    viSites1 = S0.P.miSites(:, iSite1);
else
    viSites1 = iSite1;
end
% mrSpkWav1 = tnWav2uV_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, viSites1));
mrSpkWav1 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, viSites1));
mrSpkWav1 = reshape(mrSpkWav1, [], size(mrSpkWav1,3));
vlSpkIn = auto_split_wav_(mrSpkWav1);

hFigTemp = gcf;
try close(hMsg); catch; end
if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes'), 'Yes')
    close(hFigTemp);
    split_clu_(iClu1, vlSpkIn);    
else
    close(hFigTemp);
end
end %func


%--------------------------------------------------------------------------
function edit_prm_(hObject, event)
% Edit prm file

P = get0_('P');
edit(P.vcFile_prm);
end %func


%--------------------------------------------------------------------------
function export_csv_(varargin)
% export_csv_(hObject, event)
% if nargin<2, 
fZeroIndex = 0; %zero-based index export (disable to export as matlab index starting from 1)

% S0 = get(0, 'UserData');
if nargin==2
    [S0, P, S_clu] = get0_();
elseif nargin==1
    P = varargin{1};
    vcFile_prm = P.vcFile_prm;
    S0 = load_cached_(P, 0);
    if isempty(S0), fprintf(2, 'Cannot find _jrc.mat.\n'); return; end %exit if file doesn't exist
    P = S0.P;
end

% vcFile_clu = subsFileExt(P.vcFile_prm, '_clu.mat');
% Sclu = load(vcFile_clu); %load Sclu    
% if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end

if isfield(S0, 'S_clu')
    viClu = double(S0.S_clu.viClu);
else
    fprintf(2, 'Cannot find S_clu.\n');
end
vrTime = double(S0.viTime_spk) / P.sRateHz;
viSite = double(S0.viSite_spk) - fZeroIndex; %zero base

vcFile_csv = subsFileExt_(P.vcFile_prm, '.csv');
dlmwrite(vcFile_csv, [vrTime(:), viClu(:), viSite(:)], 'precision', 9);
fprintf('wrote to %s\n', vcFile_csv);
end %func


%--------------------------------------------------------------------------
function hMsgbox = msgbox_(csMsg, fBlock, fModal)
% msgbox. Don't display if fDebug_ui is set
hMsgbox = []; 
if nargin<2, fBlock = 0; end
if nargin<3, fModal = 0; end
global fDebug_ui
if fDebug_ui==1, return; end
% if get_set_([], 'fDebug_ui', 0), return; end %don't display msgbox when debugging
% if get_set_(get0_('P'), 'fVerbose', 1)    
% fprintf('%s: %s\n', datestr(now, 'HH:MM:SS'), toString_(csMsg)); 
% end
if fBlock    
    uiwait(msgbox(csMsg, 'modal'));
else
    try
        hMsgbox = msgbox(csMsg, ifeq_(fModal, 'modal', 'non-modal'));
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function unit_annotate_(hObject, event, vcLabel)
S0 = get(0, 'UserData');
S_clu = S0.S_clu;
iClu1 = S0.iCluCopy;
if ~isfield(S_clu, 'csNote_clu'), S_clu.csNote_clu = cell(S_clu.nClu, 1); end
if nargin==3
    if isempty(vcLabel), vcLabel='';
    elseif vcLabel(1) == '='
        if ~isempty(S0.iCluPaste)
            vcLabel = sprintf('=%d', S0.iCluPaste);
        else
            msgbox_('Right-click another unit to set equal to.');
            return;
            vcLabel = '';
        end
    end
    S0.S_clu.csNote_clu{iClu1} = vcLabel;    
else
    vcNote1 = S_clu.csNote_clu{iClu1};
    if isempty(vcNote1), vcNote1=''; end
    csAns = inputdlg_(sprintf('Clu%d', iClu1), 'Annotation', 1, {vcNote1});
    if isempty(csAns), return; end
    vcLabel = csAns{1};
    S0.S_clu.csNote_clu{iClu1} = vcLabel;
end
S0 = save_log_(sprintf('annotate %d %s', iClu1, vcLabel), S0);
% set(0, 'UserData', S0);
clu_info_(S0); %update label
% update cluster
end %func


%--------------------------------------------------------------------------
function reload_prm_(hObject, event)
% Edit prm file
% 2016 07 06
[S0, P] = get0_();
S0.P = loadParam_(P.vcFile_prm);
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
function reset_position_(hObject, event)
% bottom to top left to right
S0 = get(0, 'UserData');
% P = S0.P;
for iFig=1:numel(S0.csFig)
    hFig1 = get_fig_cache_(S0.csFig{iFig});
    if ishandle(hFig1)
        set(hFig1, 'OuterPosition', S0.cvrFigPos0{iFig});
        figure(hFig1); %bring it to foreground
    end
end
end %func


%--------------------------------------------------------------------------
% function about_()
% 7/24/17 JJJ: Updated requirements and contact info
function csAbout = about_(varargin)
[vcVer, vcDate] = jrc_version_();
csAbout = { ...            
    ''; 
    sprintf('Jun Rocket Clust %s (jrc3.m)', vcVer);
    sprintf('  Last updated on %s', vcDate);
    '  Created by James Jun (jamesjun@gmail.com)';
    '';
    'Hardware Requirements';
    '  32GB ram (or 1/4 of recording size)';
    '  NVIDIA GPU (Compute Capability 3.5+: Kepler, Maxwell or Pascal)';    
    '';
    'Software Requirements';
    '  Matlab (R2014b or higher) with Toolboxes below';
    '    Parallel Processing, Image processing, Signal Processing, Statistics and Machine Learning';
    '  CUDA version supported by Matlab prallel processing toolbox (link below)';
    '    CUDA 8.0 (R2017a,b) link: https://developer.nvidia.com/cuda-downloads';
    '    CUDA 7.5 (R2016a,b) link: https://developer.nvidia.com/cuda-75-downloads-archive';
    '    CUDA 7.0 (R2015b) link: https://developer.nvidia.com/cuda-toolkit-70';
    '    CUDA 6.5 (R2015a) link: https://developer.nvidia.com/cuda-toolkit-65';
    '    CUDA 6.0 (R2014b) link: https://developer.nvidia.com/cuda-toolkit-60';
    '  Latest NVidia GPU driver (link below)';
    '    http://www.nvidia.com/Download/index.aspx';
    '  Visual Studio 2013 Express (link below)';
    '    https://www.dropbox.com/s/jrpmto1mwdp4uga/en_visual_studio_express_2013_for_windows_desktop_with_update_5_x86_web_installer_6815514.exe?dl=1'
}; 
if nargout==0, disp_cs_(csAbout); end
end


%--------------------------------------------------------------------------
function keyPressFcn_FigWavCor_(hObject, event)
S0 = get(0, 'UserData');
switch lower(event.Key)
    case 'm' %merge
        ui_merge_();
    case 's' %split
        auto_split_(1); %multi
    case {'d', 'backspace', 'delete'} %delete
        ui_delete_();        
end %switch
end %func


%--------------------------------------------------------------------------
function help_FigWav_(hObject, event)
[~, S_fig] = get_fig_cache_('FigWav');
msgbox_(S_fig.csHelp, 1);
end %func


%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
function S = select_polygon_(hPlot)
S = get(hPlot, 'UserData');
% try delete(S.hPlot_split); catch; end
update_plot2_proj_();
% try delete(S.hPoly); catch; end

S.hPoly = impoly_(); %get a polygon drawing from user
if isempty(S.hPoly), S=[]; return; end;
mrPolyPos = getPosition(S.hPoly);

vrXp = get(hPlot, 'XData');
vrYp = get(hPlot, 'YData');
vlKeep1 = inpolygon(vrXp, vrYp, mrPolyPos(:,1), mrPolyPos(:,2));
% viKeep1 = find(inpoly([vrXp(:), vrYp(:)], mrPolyPos));

[viEvtPlot,~,~] = ind2sub(S.tr_dim, S.viPlot); 
viEvtKeep1 = unique(viEvtPlot(vlKeep1));
viEvtKeep = find(ismember(viEvtPlot, viEvtKeep1));

update_plot2_proj_(vrXp(viEvtKeep), vrYp(viEvtKeep));
% S.hPlot_split = line(vrXp(viEvtKeep), vrYp(viEvtKeep), 'Color', [1 0 0], 'Marker', 'o', 'MarkerSize', 1, 'LineStyle', 'none'); 
% nSites = S.tr_dim(2);
% axis([0 nSites 0 nSites]); %reset view

set(hPlot, 'UserData', S);
hold off;
end %fund


%--------------------------------------------------------------------------
function [S_clu, iClu2] = clu_reorder_(S_clu, iClu1, iClu2)
% Move iClu2 next to iClu1
% from end to location. iClu1: place next to iClu1
% reorder clu location
% cluster is already appended
if nargin<2,  iClu2 = S_clu.nClu; end

if nargin < 2
    iClu1 = find(S_clu.viSite_clu < S_clu.viSite_clu(end), 1, 'last');
    if isempty(iClu1), return; end
end
iClu2 = iClu1+1; %move one right to iClu1
if iClu2 == S_clu.nClu, return; end %if no change in position return

vlAdd = S_clu.viClu>iClu1 & S_clu.viClu < S_clu.nClu;
S_clu.viClu(S_clu.viClu == S_clu.nClu) = iClu2;
S_clu.viClu(vlAdd) = S_clu.viClu(vlAdd) + 1;

viMap_clu = 1:S_clu.nClu;
viMap_clu(iClu2) = S_clu.nClu;
viMap_clu(iClu1+2:end) = (iClu1+1):(S_clu.nClu-1);

% S_clu = struct_reorder_(S_clu, viMap_clu, ...
%     'cviSpk_clu', 'vnSpk_clu', 'viSite_clu', 'csNote_clu', 'vrPosX_clu', 'vrPosY_clu');
S_clu = S_clu_select_(S_clu, viMap_clu);
end %func


%--------------------------------------------------------------------------
function export_mrWav_clu_(h,e)
% Export selected cluster waveforms (iCopy, iPaste) set in the GUI

S0 = get(0, 'UserData');
P = S0.P;
S_clu = S0.S_clu;
viSite1 = P.miSites(:,S_clu.viSite_clu(S0.iCluCopy));
mrWav_clu1 = S_clu.tmrWav_clu(:,viSite1,S0.iCluCopy);
eval(sprintf('mrWav_clu%d = mrWav_clu1;', S0.iCluCopy));
eval(sprintf('assignWorkspace_(mrWav_clu%d);', S0.iCluCopy));
if ~isempty(S0.iCluPaste);
    mrWav_clu2 = S_clu.tmrWav_clu(:,viSite1,S0.iCluPaste);
    eval(sprintf('mrWav_clu%d = mrWav_clu2;', S0.iCluPaste));
    eval(sprintf('assignWorkspace_(mrWav_clu%d);', S0.iCluPaste));
end
end %func


%--------------------------------------------------------------------------
function export_tmrWav_clu_(hObject, event)
% exports spike waveforms by clusters
[S_clu, P] = get0_('S_clu', 'P');
tmrWav_clu = S_clu.tmrWav_clu;
assignWorkspace_(tmrWav_clu);
end %func


%--------------------------------------------------------------------------
function export_tnWav_spk_(h,e)
% Export all spike waveforms from selected cluster

global tnWav_spk tnWav_raw
S0 = get(0, 'UserData');
P = S0.P;
S_clu = S0.S_clu;
iClu1 = S0.iCluCopy;
viSpk1 = S_clu.cviSpk_clu{iClu1};
nSpk1 = numel(viSpk1);
nSites = numel(P.viSite2Chan);
dimm_spk1 = size(tnWav_spk);    dimm_spk1(2) = nSites;  
dimm_raw1 = size(tnWav_raw);    dimm_raw1(2) = nSites;
tnWav_spk1 = zeros(dimm_spk1, 'like', tnWav_spk);
tnWav_raw1 = zeros(dimm_raw1, 'like', tnWav_raw);
miSites_spk1 = P.miSites(:, S0.viSite_spk);
for iSpk = 1:nSpk1
    tnWav_spk1(:, miSites_spk1(:,iSpk), iSpk) = tnWav_spk(:,:,iSpk);
    tnWav_raw1(:, miSites_spk1(:,iSpk), iSpk) = tnWav_raw(:,:,iSpk);
end
eval(sprintf('tnWav_spk_clu%d = tnWav_spk1;', iClu1));
eval(sprintf('assignWorkspace_(tnWav_spk_clu%d);', iClu1));
eval(sprintf('tnWav_raw_clu%d = tnWav_raw1;', iClu1));
eval(sprintf('assignWorkspace_(tnWav_raw_clu%d);', iClu1));
end %func


%--------------------------------------------------------------------------
function save_figures_(vcExt)
% bottom to top left to right
global fDebug_ui

vcPrefix = sprintf('jrc3_%s_', datestr(now, 'yymmdd-HHMM'));
csAns = inputdlg_('Figure name prefix', 'Save figure set', 1, {vcPrefix});
if isempty(csAns), return; end
vcPrefix = csAns{1};

csFig = get0_('csFig');
fprintf('Saving figures...\n'); t1=tic;
for iFig=1:numel(csFig)
%     hFig1 = P.vhFig(iField);
    hFig1 = get_fig_cache_(csFig{iFig});
    if ~ishandle(hFig1), continue; end
    vcFile1 = [vcPrefix, get(hFig1, 'Tag'), vcExt];
    if fDebug_ui==1, continue; end
%     if get_set_([], 'fDebug_ui', 0), continue; end %skip saving for debugging
    switch lower(vcExt)
        case '.fig'
            savefig(hFig1, vcFile1, 'compact');
        otherwise
            saveas(hFig1, vcFile1, vcExt(2:end));
    end
    fprintf('\t%s\n', vcFile1);
end
fprintf('\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function plot_SpikePos_(S0, event)
% global mrWav
fPlotX = 0; fPlotAllSites = 0; skip_bg = 4;
if key_modifier_(event, 'shift'), fPlotX = 1; end
if key_modifier_(event, 'alt'), fPlotAllSites = 1; end
if isempty(S0), S0 = get(0, 'UserData'); end
P = S0.P;
S_clu = S0.S_clu;
[hFig, S_fig] = get_fig_cache_('FigTime');
nSites = numel(P.viSite2Chan);
if nargin<2, fPlotX = 0; end %plot x if set

if fPlotAllSites
    viSites = 1:nSites;
else
    iSite_clu = S_clu.viSite_clu(S0.iCluCopy); %only plot spikes from site     
    viSites =  P.miSites(:, iSite_clu)';
end

% determine spike positinos
[vrX0, vrY0, vrA0, viClu0, vrT0] = get_spike_clu_(S_clu, viSites); %background
[vrX2, vrY2, vrA2, viClu2, vrT2] = deal([]);
if ~fPlotAllSites
    [vrX1, vrY1, vrA1, viClu1, vrT1] = multiindex_(find(viClu0==S0.iCluCopy), ...
        vrX0, vrY0, vrA0, viClu0, vrT0);
    if ~isempty(S0.iCluPaste)
        [vrX2, vrY2, vrA2, viClu2, vrT2] = multiindex_(find(viClu0==S0.iCluPaste), ...
            vrX0, vrY0, vrA0, viClu0, vrT0);
    end
else
    % display chain
    viClu_Chain(1) = S0.iCluCopy;
    iClu_next = get_next_clu_(S_clu, S0.iCluCopy);  
    while ~isempty(iClu_next)        
        viClu_Chain(end+1) = iClu_next;
        iClu_next = get_next_clu_(S_clu, iClu_next);        
    end
    [vrX1, vrY1, vrA1, viClu1, vrT1] = multiindex_(find(ismember(viClu0, viClu_Chain)), ...
        vrX0, vrY0, vrA0, viClu0, vrT0);
end

% if fPlot_ampDist, plot_ampDist_(cmrVpp_site, P); end

if skip_bg>1 % subsample background
    [vrA0, vrX0, vrY0, vrT0, viClu0] = multiindex_(1:skip_bg:numel(vrA0), ...
        vrA0, vrX0, vrY0, vrT0, viClu0);
end

S_fig.xylim = [get(S_fig.hAx, 'XLim'), get(S_fig.hAx, 'YLim')]; %store limits
S_fig.xylim_track = S_fig.xylim;

%------------
% Draw
if ~isfield(S_fig, 'vhAx_track')
    S_fig.vhAx_track = axes('Parent', hFig, 'Position', [.05 .2 .9 .7], 'XLimMode', 'manual', 'YLimMode', 'manual');
    hold(S_fig.vhAx_track, 'on');     
    xlabel(S_fig.vhAx_track, 'Time (s)');     
    S_fig.hPlot0_track = scatter(nan, nan, 5, nan, 'filled'); %place holder
    S_fig.hPlot1_track = line(nan, nan, 'Color', [0 0 1], 'Marker', 'o', 'MarkerSize', 5, 'LineStyle', 'none');
    S_fig.hPlot2_track = line(nan, nan, 'Color', [1 0 0], 'Marker', 'o', 'MarkerSize', 5, 'LineStyle', 'none');
else
    toggleVisible_({S_fig.vhAx_track, S_fig.hPlot0_track, S_fig.hPlot1_track, S_fig.hPlot2_track}, 1);
    set(S_fig.vhAx_track, 'Visible', 'on');    
end
% axes(S_fig.vhAx_track);
mouse_figure(hFig, S_fig.vhAx_track);
toggleVisible_({S_fig.hAx, S_fig.hRect, S_fig.hPlot0, S_fig.hPlot1, S_fig.hPlot2}, 0);
if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
toggleVisible_(S_fig.hPlot0_track, S_fig.fPlot0);
if fPlotX
    set(S_fig.hPlot0_track, 'XData', vrT0, 'YData', vrX0, 'CData', log10(vrA0));
    update_plot_(S_fig.hPlot1_track, vrT1, vrX1);
    update_plot_(S_fig.hPlot2_track, vrT2, vrX2);
    ylabel(S_fig.vhAx_track, 'X Pos [pix]');
    S_fig.xylim_track(3:4) = [0 1.5];
else
    set(S_fig.hPlot0_track, 'XData', vrT0, 'YData', vrY0, 'CData', log10(vrA0));
    update_plot_(S_fig.hPlot1_track, vrT1, vrY1);
    update_plot_(S_fig.hPlot2_track, vrT2, vrY2);
    ylabel(S_fig.vhAx_track, 'Y Pos [pix]');
    S_fig.xylim_track(3:4) = round(median(vrY1)) + [-1,1] * floor(P.maxSite);
end
axis(S_fig.vhAx_track, S_fig.xylim_track);
colormap(S_fig.vhAx_track, flipud(colormap('gray')));
set(S_fig.vhAx_track, 'CLim', [.5 3.5]);
grid(S_fig.vhAx_track, 'on');

% Set title
if isempty(S0.iCluPaste)
    vcTitle = sprintf('Color: log10 Vpp [uV]; Clu%d(black); Press [T] to return; [B]ackground', S0.iCluCopy);
else
    vcTitle = sprintf('Color: log10 Vpp [uV]; Clu%d(black); Clu%d(red); Press [T] to return; [B]ackground', S0.iCluCopy, S0.iCluPaste);
end
title(S_fig.vhAx_track, vcTitle);
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function [vrX_spk, vrY_spk, vrA_spk, viClu_spk, vrT_spk] = get_spike_clu_(S_clu, viSites, viClu_plot)
% global mrWav
% error('get_spike_clu_: resolve mrWav dependency');
if nargin<3, viClu_plot=[]; end
S0 = get(0, 'UserData');
P = S0.P;

nSpikes = numel(S0.viTime_spk);
[vrX_spk, vrY_spk, vrA_spk, viClu_spk, vrT_spk] = deal(nan(nSpikes, 1, 'single'));
for iSite = viSites
    viSite1 = P.miSites(:, iSite);
    viSpk1 = find(S0.viSite_spk == iSite);
    viClu1 = S_clu.viClu(viSpk1);
    if ~isempty(viClu_plot)        
        vl1_clu = ismember(viClu1, viClu_plot);
        viSpk1 = viSpk1(vl1_clu);
        viClu1 = viClu1(vl1_clu);
    end
    if isempty(viSpk1), continue; end
    viTime1 = S0.viTime_spk(viSpk1);    
%     [vrX_spk(viSpk1), vrY_spk(viSpk1), vrA_spk(viSpk1)] = ...
%         spikePos_(viSpk1, viSite1, P);
    viClu_spk(viSpk1) = viClu1;
    vrT_spk(viSpk1) = viTime1;
end %for

% select
vl_spk = ~isnan(vrY_spk);
vrA_spk = (vrA_spk(vl_spk));
vrY_spk = vrY_spk(vl_spk) / P.um_per_pix;
vrX_spk = vrX_spk(vl_spk) / P.um_per_pix;
vrT_spk = single(vrT_spk(vl_spk)) / P.sRateHz;
viClu_spk = viClu_spk(vl_spk);

% sort by amplitude
[vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk] = ...
    sort_ascend_(vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk);
end %func


%--------------------------------------------------------------------------
% function [vrX1_spk, vrY1_spk, vrVpp_spk, mrVpp1, trWav1] = spikePos_(viSpk1, viSite1, P)
% % Get spike position from a spike index viSpk1
% % Determine Vpp
% % viTime1 = randomSelect(viTime1, P.nShow_proj); %select fewer
% % [trWav1, mrVpp1] = mr2tr_spk_(mrWav, viTime1, viSite1, P);
% 
% fCentroid = 0;
% 
% if fCentroid
%     P.fInterp_mean = 0;
% %     [mrWav_mean1, imax_site1, trWav_int1] = interp_align_mean_(trWav1, P);
%     mrXY_spk = centroid_pca_(trWav1, P.mrSiteXY(viSite1, :));
% end
% vrX1_site = P.mrSiteXY(viSite1, 1);
% vrY1_site = P.mrSiteXY(viSite1, 2);    
% mrVpp1_sq = mrVpp1.^2;
% vrVpp1_sq_sum = sum(mrVpp1_sq);
% if ~fCentroid
%     vrX1_spk = sum(bsxfun(@times, mrVpp1_sq, vrX1_site)) ./ vrVpp1_sq_sum;
%     vrY1_spk = sum(bsxfun(@times, mrVpp1_sq, vrY1_site)) ./ vrVpp1_sq_sum;
% else
%     vrX1_spk = mrXY_spk(:,1);
%     vrY1_spk = mrXY_spk(:,1);
% end
% vrVpp_spk = sqrt(vrVpp1_sq_sum);
% end %func


%--------------------------------------------------------------------------
function update_FigTime_()
% display features in a new site

[hFig, S_fig] = get_fig_cache_('FigTime');
S0 = get(0, 'UserData');
P = S0.P;
if ~isVisible_(S_fig.hAx), return ;end
% P.vcFet_show = S_fig.csFet{S_fig.iFet};
set0_(P);
[vrFet0, vrTime0, vcYlabel] = getFet_site_(S_fig.iSite);
if ~isfield(S_fig, 'fPlot0'), S_fig.fPlot0 = 1; end
toggleVisible_(S_fig.hPlot0, S_fig.fPlot0);
update_plot_(S_fig.hPlot0, vrTime0, vrFet0);
set(S_fig.hPlot1, 'YData', getFet_site_(S_fig.iSite, S0.iCluCopy));
% set(S_fig.hPlot0, 'XData', vrTime0, 'YData', vrFet0);
if ~isempty(S0.iCluPaste)
    set(S_fig.hPlot2, 'YData', getFet_site_(S_fig.iSite, S0.iCluPaste));
else
    hide_plot_(S_fig.hPlot2);
%     set(S_fig.hPlot2, 'XData', nan, 'YData', nan);
end
% switch lower(P.vcFet_show)
%     case 'vpp'
ylim(S_fig.hAx, [0, 1] * S_fig.maxAmp);
imrect_set_(S_fig.hRect, [], [0, 1] * S_fig.maxAmp);    
%     otherwise
%         ylim(S_fig.hAx, [0, 1] * P.maxAmp);
%         imrect_set_(S_fig.hRect, [], [0, 1] * P.maxAmp);    
% end
grid(S_fig.hAx, 'on');
ylabel(S_fig.hAx, vcYlabel);
end %func


%--------------------------------------------------------------------------
function S_clu = split_clu_(iClu1, vlIn)
% split cluster. 
% global mrWav
[P, S_clu, viSite_spk] = get0_('P', 'S_clu', 'viSite_spk');
hMsg = msgbox_open_('Splitting...');
figure(get_fig_cache_('FigWav'));

% create a new cluster (add at the end)
n2 = sum(vlIn); %number of clusters to split off
iClu2 = max(S_clu.viClu) + 1;

% update cluster count and index
S_clu.nClu = double(iClu2);
S_clu.vnSpk_clu(iClu1) = S_clu.vnSpk_clu(iClu1) - n2;
S_clu.vnSpk_clu(iClu2) = sum(vlIn);
viSpk1 = find(S_clu.viClu==iClu1);
viSpk2 = viSpk1(vlIn);
viSpk1 = viSpk1(~vlIn);
S_clu.cviSpk_clu{iClu1} = viSpk1;
S_clu.cviSpk_clu{iClu2} = viSpk2;
S_clu.viSite_clu(iClu1) = mode(viSite_spk(viSpk1));
S_clu.viSite_clu(iClu2) = mode(viSite_spk(viSpk2));
try % erase annotaiton    
    S_clu.csNote_clu{iClu1} = ''; 
    S_clu.csNote_clu{end+1} = '';  %add another entry
catch
end
S_clu.viClu(viSpk2) = iClu2; %change cluster number
S_clu = S_clu_update_(S_clu, iClu1, P);
S_clu = S_clu_update_(S_clu, iClu2, P);
% S_clu = S_clu_position_(S_clu, [iClu1, iClu2]);

% Bring the new cluster right next to the old one using index swap
[S_clu, iClu2] = clu_reorder_(S_clu, iClu1); 
for iClu3 = iClu2+1:S_clu.nClu % update cluster chain info
    S_clu = S_clu_update_note_(S_clu, iClu3, get_next_clu_(S_clu, iClu3) + 1);
end

% update all the other views
[S_clu, S0] = S_clu_commit_(S_clu, 'split_clu_');
S0 = save_log_(sprintf('split %d', iClu1), S0); %@TODO: specify which cut to use
plot_FigWav_(S0); %redraw plot
plot_FigWavCor_(S0);
% select two clusters being split
button_CluWav_simulate_(iClu1, iClu2);
close_(hMsg);
fprintf('%s [W] splited Clu %d\n', datestr(now, 'HH:MM:SS'), iClu1);
end %func


%--------------------------------------------------------------------------
function [fSplit, vlIn] = plot_split_(S1)
% global tnWav_spk mrFet
% find site
% global mrWav
mrPolyPos = getPosition(S1.hPoly);
site12_show = floor(mean(mrPolyPos));
site12 = S1.viSites_show(site12_show+1);  

% get amp
S0 = get0_;
S_clu = S0.S_clu;
P = S0.P;
iClu1 = S0.iCluCopy;
% trWav12 = mr2tr_spk_(mrWav, viTime1, site12, P);
% trWav12 = bit2uV_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, site12));
% trWav12 = tnWav2uV_(tnWav_sites_(tnWav_spk, S_clu.cviSpk_clu{iClu1}, site12));
trWav12 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, site12));
if diff(site12) == 0, trWav12(:,2,:) = trWav12(:,1,:); end
vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;
switch lower(P.vcFet_show)   
    case {'vpp', 'vmin', 'vmax'}     
        mrAmin12 = abs(squeeze_(min(trWav12)))';
        mrAmax12 = abs(squeeze_(max(trWav12)))';
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (min amp)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (min amp)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (max amp)', site12(1));
            vxPoly = vxPoly; %max amp are scaled half
        end  
        
    case {'cov', 'spacetime'}   
        [mrAmin12, mrAmax12] = calc_cov_spk_(S_clu.cviSpk_clu{iClu1}, site12);
        [mrAmin12, mrAmax12] = multifun_(@(x)abs(x'), mrAmin12, mrAmax12);
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (cov1)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (cov1)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (cov2)', site12(1));
        end  
        
    case 'pca'
        [mrAmin12, mrAmax12] = pca_pc_spk_(S_clu.cviSpk_clu{iClu1}, site12);
        [mrAmin12, mrAmax12] = multifun_(@(x)abs(x'), mrAmin12, mrAmax12);
        vyPlot = mrAmin12(:,2);
        vcYlabel = sprintf('Site %d (PC1)', site12(2));
        if site12(2) > site12(1)
            vxPlot = mrAmin12(:,1);
            vcXlabel = sprintf('Site %d (PC1)', site12(1));
        else
            vxPlot = mrAmax12(:,1);
            vcXlabel = sprintf('Site %d (PC2)', site12(1));
        end 
        
    otherwise 
        error('plot_split: vcFetShow: not implemented');
%         vxPoly = (mrPolyPos([1:end,1],1) - site12_show(1)) * S1.maxAmp;
%         vyPoly = (mrPolyPos([1:end,1],2) - site12_show(2)) * S1.maxAmp;        
%         trFet12 = trFet_([], site12, S_clu.cviSpk_clu{iClu1});
%         vyPlot = squeeze_(trFet12(1, 2, :));
%         vcYlabel = sprintf('Site %d (%s1)', site12(2), P.vcFet);
%         if site12(2) > site12(1)
%             vxPlot = squeeze_(trFet12(1, 1, :));
%             vcXlabel = sprintf('Site %d (%s1)', site12(1), P.vcFet);
%         else
%             vxPlot = squeeze_(trFet12(min(2,P.nPcPerChan), 1, :));
%             vcXlabel = sprintf('Site %d (%s2)', site12(1), P.vcFet);
%         end       
end

vlIn = inpolygon(vxPlot, vyPlot, vxPoly, vyPoly);

% Plot temporary figure (auto-close)
hFig = figure(10221); clf;
resize_figure_(hFig, [0 0 .5 1]);
subplot(2,2,[1,3]); hold on;
line(vxPlot, vyPlot, 'Color', P.mrColor_proj(2,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
hPlot = line(vxPlot(vlIn), vyPlot(vlIn), 'Color', P.mrColor_proj(3,:), 'Marker', 'o', 'MarkerSize', 2, 'LineStyle', 'none');
% plot(vxPoly, vyPoly, 'b+-'); %boundary
title(sprintf('Cluster %d (%d spikes)', iClu1, S_clu.vnSpk_clu(iClu1)));
xlabel(sprintf('Site %d', site12(1)));
ylabel(vcYlabel);   xlabel(vcXlabel);
grid on;

% user edit polygon
hPoly = impoly_(gca, [vxPoly(:), vyPoly(:)]);
hFunc = makeConstrainToRectFcn('impoly',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(hPoly, hFunc); 

hMsgbox = msgbox_('Press OK after adjusting polygon', 0);
uiwait(hMsgbox);

vlIn = poly_mask_(hPoly, vxPlot, vyPlot);
set(hPlot, 'XData', vxPlot(vlIn), 'YData',vyPlot(vlIn));
% if P.fWav_raw_show
%     trWav12_raw = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu1}, site12, 1));
%     mrWavX = squeeze_(trWav12_raw(:, 1, :));
%     mrWavY = squeeze_(trWav12_raw(:, 2, :));
% else
    mrWavX = squeeze_(trWav12(:, 1, :));
    mrWavY = squeeze_(trWav12(:, 2, :));
% end        
vrT = (P.spkLim(1):P.spkLim(end)) / P.sRateHz * 1000;
viIn = randomSelect_(find(vlIn), P.nSpk_show);
viOut = randomSelect_(find(~vlIn), P.nSpk_show);

subplot 222; hold on;
plot(vrT, mrWavX(:,viOut), 'k', vrT, mrWavX(:,viIn), 'r');
title(vcXlabel); ylabel('Voltage (\muV)'); xlabel('Time (ms)');
grid on;

subplot 224; hold on;
plot(vrT, mrWavY(:,viOut), 'k', vrT, mrWavY(:,viIn), 'r');
title(vcYlabel); ylabel('Voltage (\muV)'); xlabel('Time (ms)');
grid on;

if strcmpi(questdlg_('Split?', 'Confirmation', 'Yes', 'No', 'Yes'), 'Yes')    
    fSplit = 1;
else
    fSplit = 0;
end
try close(hFig); catch; end
end %func


%--------------------------------------------------------------------------
function update_spikes_(varargin)
S0 = get(0, 'UserData');
hMsg = msgbox_open_('Updating spikes');
fig_prev = gcf;
figure_wait_(1);
[~, S_fig] = get_fig_cache_('FigWav');
% plot_tnWav_clu_(S_fig, S0.P); %do this after plotSpk_
plot_spkwav_(S_fig, S0.P);
close_(hMsg);
figure_wait_(0);
figure(fig_prev);
end %func


%--------------------------------------------------------------------------
% function tnWav1 = tnWav_clu_site_(iClu1, viSite)
% % time x site x spikes
% % Return waveforms for cluster and specific site #
% global tnWav_spk tnWav_raw
% if nargin<2, nSpk_show=inf; end
% S0 = get(0, 'UserData');
% P = S0.P;
% [viSpk_clu1, viiSpk_clu1] = S_clu_viSpk_(iClu1, S0);
% 
% % trim
% viSpk_clu1 = randomSelect_(viSpk_clu1, P.nSpk_show);
% 
% if P.fWav_raw_show
%     trWav1 = raw2uV_(tnWav_raw(:,:,viSpk_clu1), P);
% else
%     trWav1 = bit2uV_(tnWav_spk(:,:,viSpk_clu1), P);
% end
% end %func


%--------------------------------------------------------------------------
function [mrPos_spk, mrVpp_spk] = tnWav_centroid_(tnWav_spk, viSite_spk, P)
% 2 x nSpk matrix containing centroid and amplitude

mrVpp_spk = trWav2fet_(tnWav_spk, P); % apply car
% mrVpp_spk = tr2Vpp(tnWav_spk, P)';

mrVpp_spk1 = mrVpp_spk .^ 2; % compute covariance with the center site
mrVpp_spk1 = bsxfun(@rdivide, mrVpp_spk1, sum(mrVpp_spk1, 1));
miSite_spk = P.miSites(:, viSite_spk);
mrSiteXY = single(P.mrSiteXY);
mrSiteX_spk = reshape(mrSiteXY(miSite_spk(:), 1), size(miSite_spk));
mrSiteY_spk = reshape(mrSiteXY(miSite_spk(:), 2), size(miSite_spk));
mrPos_spk = [sum(mrVpp_spk1 .* mrSiteX_spk, 1); sum(mrVpp_spk1 .* mrSiteY_spk, 1)];   
end %func


%--------------------------------------------------------------------------
function [mrPos_spk, viSpk_re] = position_spk_(viSite_spk, tnWav_spk, P)

mrPos_site1 = P.mrSiteXY(P.miSites(1, viSite_spk), :)'; %first pos

% determine centroid location and second largest amplitude
[mrPos_spk, mrA_spk] = tnWav_centroid_(tnWav_spk, viSite_spk, P);
[~, viiSite2_spk] = max(mrA_spk((2:end), :)); %find second max
miSites2 = P.miSites(2:end, :);
viiSite2_spk = sub2ind(size(miSites2), viiSite2_spk(:), viSite_spk(:));
viSite_spk2 = miSites2(viiSite2_spk);

% Find where second largest site is closer to the spike centroid
mrPos_site2 = P.mrSiteXY(viSite_spk2, :)';
dist__ = @(mr1,mr2)sum((mr1-mr2).^2);
viSpk_re = find(dist__(mrPos_spk,mrPos_site2) < dist__(mrPos_spk, mrPos_site1));
end %func


%--------------------------------------------------------------------------
function viTime1 = recenter_spk_(mrWav, viTime, viSite, P)
spkLim = [-1,1] * abs(P.spkLim(1));
viTime0 = [spkLim(1):spkLim(end)]'; %column
miTime = bsxfun(@plus, int32(viTime0), int32(viTime(:)'));
miTime = min(max(miTime, 1), size(mrWav, 1));
miSite = repmat(viSite(:)', numel(viTime0), 1); 
mrWav_spk = mrWav(sub2ind(size(mrWav), miTime, miSite));
[~, viMin_spk] = min(mrWav_spk,[],1);
viTime_off = int32(gather_(viMin_spk') + spkLim(1) - 1);
viTime1 = viTime + viTime_off;
% disp(mean(viTime_off~=0))
end %func


%--------------------------------------------------------------------------
function hPatch = plot_probe_(mrSiteXY, vrSiteHW, viSite2Chan, vrVpp, hFig)
if nargin<3, viSite2Chan=[]; end
if nargin<4, vrVpp=[]; end
if nargin<5, hFig=[]; end

if isempty(hFig)
    hFig = gcf;
else
    figure(hFig);
end

vrX = [0 0 1 1] * vrSiteHW(2); 
vrY = [0 1 1 0] * vrSiteHW(1);

mrPatchX = bsxfun(@plus, mrSiteXY(:,1)', vrX(:));
mrPatchY = bsxfun(@plus, mrSiteXY(:,2)', vrY(:));
nSites = size(mrSiteXY,1);    
if ~isempty(vrVpp)
    hPatch = patch(mrPatchX, mrPatchY, repmat(vrVpp(:)', [4, 1]), 'EdgeColor', 'k', 'FaceColor', 'flat'); %[0 0 0], 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', [0 0 0], 'FaceAlpha', 0); 
    caxis([0, max(vrVpp)]);
    colormap jet;
else
    hPatch = patch(mrPatchX, mrPatchY, 'w', 'EdgeColor', 'k'); %[0 0 0], 'EdgeColor', 'none', 'FaceColor', 'flat', 'FaceVertexCData', [0 0 0], 'FaceAlpha', 0); 
end
if ~isempty(viSite2Chan)
    csText = arrayfun(@(i)sprintf('%d/%d', i, viSite2Chan(i)), 1:numel(viSite2Chan), 'UniformOutput', 0);
else
    csText = arrayfun(@(i)sprintf('%d', i), 1:nSites, 'UniformOutput', 0);
end
hText = text(mrSiteXY(:,1), mrSiteXY(:,2), csText, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left');
axis([min(mrPatchX(:)), max(mrPatchX(:)), min(mrPatchY(:)), max(mrPatchY(:))]);
title('Site# / Chan# (zoom: wheel; pan: hold wheel & drag)');
% vrPos = get(gcf, 'Position');
xlabel('X Position (\mum)');
ylabel('Y Position (\mum)');
mouse_figure(hFig);
end


%--------------------------------------------------------------------------
function S_clu = post_merge_shank_(cS_clu_shank, P)
% combine cluster info divided by clusters
% sort clusters according to their position (x,y combination)
nRepeat_merge = 4;
nShanks = numel(cS_clu_shank);
nClu = 0;
t_runtime = 0;
S0 = get(0, 'UserData');
nSpk = numel(S0.viTime_spk);
nTime_clu = get_(P, 'nTime_clu');
for iShank = 1:nShanks
    S_clu1 = postCluster_(cS_clu_shank{iShank}, P);
%     if nTime_clu > 1
%         S_clu1 = post_merge_wav_(S_clu_refresh_(S_clu1), nRepeat_merge);    
%     end
    cS_clu_shank{iShank} = S_clu1;        
    if iShank==1, viClu = zeros(nSpk, 1, 'like', S_clu1.viClu); end
    vl_zero1 = find(S_clu1.viClu>0);
    viClu(S_clu1.viSpk_shank(vl_zero1)) = S_clu1.viClu(vl_zero1) + nClu;
    
    % accum
    nClu = nClu + S_clu1.nClu;
    t_runtime = t_runtime + S_clu1.t_runtime;
end %shank

% output
S_clu = makeStruct_(nClu, cS_clu_shank, P, t_runtime, viClu); % add up number of clu
end %func


%--------------------------------------------------------------------------
function S_clu = postCluster_(S_clu, P)
% fSortClu = 0; %debug
if isfield(S_clu, 'viClu')
    S_clu = rmfield(S_clu, 'viClu');
end
if isempty(S_clu), return; end

switch lower(P.vcDetrend_postclu)
%         case {'hidehiko', 'hide'}
%             S_clu.icl = selec_rho_delta_with_slope(S_clu, P.delta1_cut);
    case 'local' % 
        S_clu.icl = detrend_local_(S_clu, P, 1);
    case 'none'
        S_clu.icl = find(S_clu.rho(:) > 10^(P.rho_cut) & S_clu.delta(:) > 10^(P.delta1_cut));
    case 'global'
%         S_clu.icl = detrend_ztran_(S_clu.rho, S_clu.delta, P.rho_cut, P.delta1_cut);
        S_clu.icl = detrend_local_(S_clu, P, 0);
    case 'logz'
        S_clu.icl = log_ztran_(S_clu.rho, S_clu.delta, P.rho_cut, 4+P.delta1_cut);
    otherwise
        fprintf(2, 'postCluster_: vcDetrend_postclu = ''%s''; not supported.\n', P.vcDetrend_postclu);
end
% end

% Update P
S_clu.P.min_count = P.min_count;
S_clu.P.delta1_cut = P.delta1_cut;
S_clu.P.rho_cut = P.rho_cut;
S_clu.viClu = [];
S_clu = assign_clu_count_(S_clu, P); % enforce min count algorithm    

% debug output
if nargout==0
    vrXp = log10(S_clu.rho);
    vrYp = log10(S_clu.delta);
    figure; hold on; 
    plot(vrXp, vrYp, 'b.', vrXp(S_clu.icl), vrYp(S_clu.icl), 'ro');
    plot(get(gca, 'XLim'), P.delta1_cut*[1 1], 'r-');
    plot(P.rho_cut*[1 1], get(gca, 'YLim'), 'r-');
    xlabel('log10 rho');
    ylabel('log10 delta');    
end 
end %func


%--------------------------------------------------------------------------
function [icl, x, y] = log_ztran_(x, y, x_cut, y_cut)
% [icl, x, y] = detrend_ztran_(x, y, x_cut, y_cut)
% [icl, x, y] = detrend_ztran_(x, y, n_icl)
if nargin == 3
    n_icl = x_cut;
else
    n_icl = [];
end

x = log10(x(:));
y = log10(y(:));

vlValid = isfinite(x) & isfinite(y);
y(vlValid) = zscore_(y(vlValid));

if isempty(n_icl)
    icl = find(x>=x_cut & y>=y_cut);
else
    [~, ix] = sort(y(vlValid), 'descend');   
    viValid = find(vlValid);
    icl = viValid(ix(1:n_icl));
end

if nargout==0
    figure; plot(x,y,'.', x(icl),y(icl),'ro'); grid on; 
end
end %func


%--------------------------------------------------------------------------
function [icl, x, z] = detrend_ztran_(x0, y0, x_cut, y_cut)
% [icl, x, y] = detrend_ztran_(x, y, x_cut, y_cut)
% [icl, x, y] = detrend_ztran_(x, y, n_icl)
if nargin == 3
    n_icl = x_cut;
else
    n_icl = [];
end
y_thresh_mad = 20;
x = log10(x0(:));
y = (y0(:));
vlValid_x = x>=x_cut & x <= -1;
z_max = 100;

% compute y_mad
y_med = median(y(vlValid_x));
y_mad = y - y_med;
y_mad = y_mad / median(y_mad);
vlValid_y = abs(y_mad) < y_thresh_mad;
% vl_y0 = y0<=0;
% y = y0;
% y(vl_y0) = min(y(~vl_y0));
% y = log10(y0(:));
% max_y = max(y);

vlValid = isfinite(x) & isfinite(y);
% viDetrend = find(vlValid  & (y < max_y/4) & (x >= x_cut));
viDetrend = find(vlValid  & vlValid_y & vlValid_x);
x1 = x(viDetrend);
y1 = y(viDetrend);
xx1 = [x1+eps('single'), ones(numel(x1),1)];
m = xx1 \ y1; %determine slope based on acceptable units

xx = [x+eps('single'), ones(numel(x),1)];
y2 = y - xx * m; %detrend data
mu2 = mean(y2(viDetrend));
sd2 = std(y2(viDetrend));
z = (y2 - mu2) / sd2; %z transformation
z(z>z_max) = z_max;

if isempty(n_icl)
    icl = find(x>=x_cut & z>=10^y_cut);
else
    [~, ix] = sort(z(vlValid), 'descend');   
    viValid = find(vlValid);
    icl = viValid(ix(1:n_icl));
end

if nargout==0
    figure; plot(x,z,'.', x(icl),z(icl),'ro'); grid on; 
    title(sprintf('%d clu', numel(icl)));
end
end %func


%--------------------------------------------------------------------------
function [cl, icl] = assignCluster_(cl, ordrho, nneigh, icl)
ND = numel(ordrho);
nClu = numel(icl);

% global cluster assignment
% fprintf('Assigning clusters: '); t1 = tic;
if isempty(cl)
    cl = zeros([ND, 1], 'int32');
    cl(icl) = 1:nClu;
end

if numel(icl) == 0 || numel(icl) == 1
    cl = ones([ND, 1], 'int32');
    icl = ordrho(1);
else
    nneigh1 = nneigh(ordrho);
    for i=1:10
        vi = find(cl(ordrho)<=0);
        if isempty(vi), break; end
        vi=vi(:)';        
        for ii = vi
            cl(ordrho(ii)) = cl(nneigh1(ii));        
        end
        n1 = sum(cl<=0);
        if n1==0, break; end
        fprintf('i:%d, n0=%d, ', i, n1);
    end
    cl(cl<=0) = 1; %background
end
end %func


%--------------------------------------------------------------------------
function S_clu = assign_clu_count_(S_clu, P)
nRepeat_max = 1000;
if isempty(P.min_count), P.min_count = 0; end
if ~isfield(S_clu, 'viClu'), S_clu.viClu = []; end
if isempty(S_clu.viClu)
    nClu_pre = [];
else
    nClu_pre = S_clu.nClu;
end
nClu_rm = 0;
fprintf('assigning clusters, nClu:%d\n', numel(S_clu.icl)); t1=tic;

% fReassign = 0;
min_rho = -inf;
for iRepeat=1:nRepeat_max % repeat 1000 times max  
    S_clu.icl(S_clu.rho(S_clu.icl) < min_rho) = [];
    [S_clu.viClu, S_clu.icl] = assignCluster_(S_clu.viClu, S_clu.ordrho, S_clu.nneigh, S_clu.icl);   
    S_clu.viClu(S_clu.rho < min_rho) = 0; %noise assignment
    if isempty(P.min_count), P.min_count = 0; end
    P.min_count = max(P.min_count, S_clu.trFet_dim(1)*2);
    % http://scikit-learn.org/stable/modules/lda_qda.html

    S_clu = S_clu_refresh_(S_clu);
    
    % remove clusters unused
    viCluKill = find(S_clu.vnSpk_clu <= P.min_count);
    if isempty(viCluKill), break; end
    S_clu.icl(viCluKill) = []; 
    S_clu.viClu=[];
    nClu_rm = nClu_rm + numel(viCluKill);
    if iRepeat==nRepeat_max
        fprintf(2, 'assign_clu_count_: exceeded nRepeat_max=%d\n', nRepeat_max);
    end
end

fprintf('\n\ttook %0.1fs. Removed %d clusters having <%d spikes: %d->%d\n', ...
    toc(t1), nClu_rm, P.min_count, nClu_pre, S_clu.nClu);
end


%--------------------------------------------------------------------------
function vlIn = poly_mask_(hPoly, vxPlot, vyPlot)
mrPolyPos = getPosition(hPoly);
vxPoly = mrPolyPos([1:end,1],1);
vyPoly = mrPolyPos([1:end,1],2);
vlIn = inpolygon(vxPlot, vyPlot, vxPoly, vyPoly);
end %func


%--------------------------------------------------------------------------
function varargout = sort_ascend_(varargin)
% sort all the other fields basedon the first field in ascending order
% [a', b', c'] = sort_ascend_(a, b, c)

[varargout{1}, viSrt] = sort(varargin{1}, 'ascend');
for i=2:nargin
    varargout{i} = varargin{i}(viSrt);
end
end %func


%--------------------------------------------------------------------------
% function mrXY_spk = centroid_pca_(trWav1, mrSiteXY1, mrPv1)
% % mrXY_spk = centroid_pca_(iClu) % return 
% % mrXY_spk = centroid_pca_(trWav, mrSiteXY, mrPv)
% % mrXY_spk = centroid_pca_(mrWav, mrSiteXY, mrPv)
% P = get0_('P');
% if nargin == 1
%     fInterp_align = 0;
%     iClu1 = trWav1;
%     [~, ~, ~, mrSiteXY1, trWav1] = mrWav_int_mean_clu_(iClu1, fInterp_align);
%     mrSiteXY1 = mrSiteXY1 / P.um_per_pix; %pix unit
% end
% if ismatrix(trWav1)
%     mrWav1 = trWav1;    
% else
%     mrWav1 = mean(trWav1, 3);
% end
% if nargin<3
%     [~, mrPv1, ~] = pca(mrWav1, 'NumComponents', P.nPc_dip, 'Centered', 0); 
% end
% if ~ismatrix(trWav1)
%     mrWav1 = reshape(trWav1, size(trWav1,1), []);
% end
% mrPc_spk_sq1 = ((norm_mr_(mrPv1(:,1)) \ mrWav1)') .^ 2;
% if ~isempty(trWav1)
%     mrPc_spk_sq1 = reshape(mrPc_spk_sq1, size(trWav1,2), [])';
% end
% mrXY_spk = [sum(bsxfun(@times, mrPc_spk_sq1, mrSiteXY1(:,2)'), 2), ...
%             sum(bsxfun(@times, mrPc_spk_sq1, mrSiteXY1(:,1)'), 2)];
% mrXY_spk = bsxfun(@rdivide, mrXY_spk, sum(mrPc_spk_sq1, 2));
% end %func


%--------------------------------------------------------------------------
function [mr, vr] = norm_mr_(mr)
% Normalize each columns
vr = sqrt(sum(mr .^ 2));
mr = bsxfun(@rdivide, mr, vr);
end %func


%--------------------------------------------------------------------------
% function [mrWav_mean_int1, viSite1, sd_wav, mrSiteXY1, trWav1] = mrWav_int_mean_clu_(iClu1, fInterp_align)
% error('mrWav_int_mean_clu_: remove dependency to mrWav'); 
% global mrWav
% 
% [S0, P, S_clu] = get0_();
% if nargin<2, fInterp_align=0; end
% [miSites1, ~] = findNearSites_(P.mrSiteXY, P.maxSite_dip);  
% 
% viSpk1 = S_clu.cviSpk_clu{iClu1};
% iSite0 = S_clu.viSite_clu(iClu1); %center site        
% viSite1 = miSites1(:, iSite0);
% trWav1 = mr2tr_spk_(mrWav, S_clu.viTime(viSpk1), viSite1, P);   
% P.fInterp_align = fInterp_align; 
% [mrWav_mean_int1, imax_site_int1] = interp_align_mean_(trWav1, P);
% sd_wav = mean(std(trWav1,1,3));
% mrSiteXY1 = P.mrSiteXY(viSite1,:);
% end %func


%--------------------------------------------------------------------------
function S_clu = clu2wav_(S_clu, viSite_spk, tnWav_spk, tnWav_raw, viClu_update)
% average cluster waveforms and determine the center
% only use the centered spikes 
% viClu_copy: waveforms not changing

if nargin<4, tnWav_raw = []; end
if nargin<5, viClu_update = []; end
P = S_clu.P;
P.fMeanSubt = 0;
fVerbose = isempty(viClu_update);

if fVerbose, fprintf('Calculating cluster mean waveform.\n\t'); t1 = tic; end
if isfield(S_clu, 'nClu')
    nClu = S_clu.nClu;
else
    nClu = max(S_clu.viClu);
end
nSamples = size(tnWav_spk,1);
nSites = numel(P.viSite2Chan);
nSites_spk = size(tnWav_spk,2); % n sites per event group (maxSite*2+1);

% Prepare cluster loop
trWav_spk_clu = zeros([nSamples, nSites_spk, nClu], 'single');
vrFracCenter_clu = zeros(nClu, 1);
% miSites_clu = S_clu.P.miSites(:,S_clu.viSite_clu);
tmrWav_spk_clu = zeros(nSamples, nSites, nClu, 'single');
if ~isempty(tnWav_raw)
    nSamples_raw = size(tnWav_raw,1);
    trWav_raw_clu = zeros(nSamples_raw, nSites_spk, 'single'); 
    tmrWav_raw_clu = zeros(nSamples_raw, nSites, nClu, 'single');
else
    trWav_raw_clu = [];
    tmrWav_raw_clu = [];
end
if ~isempty(viClu_update)   
    vlClu_update = false(nClu, 1);
    vlClu_update(viClu_update) = 1;
    nClu_pre = size(S_clu.trWav_spk_clu, 3);
    vlClu_update((1:nClu) > nClu_pre) = 1;
else
    vlClu_update = true(nClu, 1);
end
for iClu=1:nClu       
    if vlClu_update(iClu)
        [viSpk1_clu1, viSite_clu1, vlSpk_clu1] = S_clu_subsample_spk_(S_clu, iClu);
        if isempty(viSpk1_clu1), continue; end        
        vrFracCenter_clu(iClu) = mean(vlSpk_clu1);        
        trWav_spk_clu(:,:,iClu) = bit2uV_(mean(single(tnWav_spk(:,:,viSpk1_clu1)),3), P);
        if P.fMeanSubt, trWav_spk_clu(:,:,iClu) = meanSubt_(trWav_spk_clu(:,:,iClu)); end
        tmrWav_spk_clu(:,viSite_clu1,iClu) = trWav_spk_clu(:,:,iClu);
        if ~isempty(tnWav_raw)
            trWav_raw_clu(:,:,iClu) = mean_tnWav_raw_(tnWav_raw(:,:,viSpk1_clu1), P);
            tmrWav_raw_clu(:,viSite_clu1,iClu) = trWav_raw_clu(:,:,iClu);
        end
    else
        trWav_spk_clu(:,:,iClu) = S_clu.trWav_spk_clu(:,:,iClu);
        tmrWav_spk_clu(:,:,iClu) = S_clu.tmrWav_spk_clu(:,:,iClu);
        if ~isempty(tnWav_raw)
            trWav_raw_clu(:,:,iClu) = S_clu.trWav_raw_clu(:,:,iClu);
            tmrWav_raw_clu(:,:,iClu) = S_clu.tmrWav_raw_clu(:,:,iClu);
        end  
    end
    if fVerbose, fprintf('.'); end
end %clu
% tmrWav_clu = meanSubt_(cumsum(tmrWav_spk_clu)); %meanSubt_ after or before?
tmrWav_clu = tmrWav_spk_clu; %meanSubt_ after or before?

% measure waveforms
[vrVmin_clu, viSite_min_clu] = min(permute(min(trWav_spk_clu),[2,3,1]),[],1);
vrVmin_clu = abs(vrVmin_clu(:));
viSite_min_clu = viSite_min_clu(:);

S_clu = struct_add_(S_clu, vrVmin_clu, viSite_min_clu, ...
    trWav_spk_clu, tmrWav_spk_clu, trWav_raw_clu, tmrWav_raw_clu, tmrWav_clu);
if fVerbose, fprintf('\n\ttook %0.1fs\n', toc(t1)); end
end %func


%--------------------------------------------------------------------------
function [viSpk_clu1, viSite_clu1, vlSpk_clu1] = S_clu_subsample_spk_(S_clu, iClu)
% subsample spikes from the requested cluster centered at the center site and mid-time range (drift)

fSelect_mid = 1;
MAX_SAMPLES = 1000;

[P, viSite_spk] = get0_('P', 'viSite_spk');
[viSpk_clu1, viSite_clu1, vlSpk_clu1] = deal([]);
% Subselect based on the center site
viSpk_clu1 = S_clu.cviSpk_clu{iClu};% 
if isempty(viSpk_clu1), return; end
iSite_clu1 = S_clu.viSite_clu(iClu);        
vlSpk_clu1 = iSite_clu1 == viSite_spk(viSpk_clu1);
viSite_clu1 = S_clu.P.miSites(:,iSite_clu1);
viSpk_clu1 = viSpk_clu1(vlSpk_clu1);
if isempty(viSpk_clu1), return; end

if fSelect_mid  
    viSpk_clu1 = spk_select_mid_(viSpk_clu1, P);
end
viSpk_clu1 = subsample_vr_(viSpk_clu1, MAX_SAMPLES);
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_sort_(S_clu, vcField_sort)
% sort clusters by the centroid position
% vcField_sort: {'', 'vrPosY_clu + vrPosX_clu'}

if nargin<2, vcField_sort = ''; end

% Sort clusters by its sites
if isempty(vcField_sort), vcField_sort = 'viSite_clu'; end

switch vcField_sort
    case 'vrPosY_clu + vrPosX_clu'
        [~, viCluSort] = sort(S_clu.vrPosY_clu + S_clu.vrPosX_clu, 'ascend');
    otherwise
        [~, viCluSort] = sort(S_clu.(vcField_sort), 'ascend');
end
S_clu.viClu = mapIndex_(S_clu.viClu, viCluSort);
S_clu = struct_reorder_(S_clu, viCluSort, ...
    'cviSpk_clu', 'vrPosX_clu', 'vrPosY_clu', 'vnSpk_clu', 'viSite_clu', 'cviTime_clu', 'csNote_clu');
% if isfield(S_clu, 'tmrWav_clu')
%     S_clu.tmrWav_clu = S_clu.tmrWav_clu(:, :, viCluSort);
% end
S_clu = S_clu_refresh_(S_clu);
end %func


%--------------------------------------------------------------------------
function [S_clu, nRemoved] = S_clu_refrac_(S_clu, P, iClu1)
% clu_refrac(Sclu, P)   %process refrac on all clusters
% clu_refrac(Sclu, P, iClu1) %process on specific clusters
% P.nSkip_refrac = 4; 
% P.fShow_refrac = 0;
viTime_spk = get0_('viTime_spk');
% remove refractory spikes
if nargin==2
%     P = varargin{1}; %second input
    nClu = max(S_clu.viClu);
    P.fShow_refrac = 1;
    nRemoved = 0;
    for iClu=1:nClu
        [S_clu, nRemoved1] = S_clu_refrac_(S_clu, P, iClu);
        nRemoved = nRemoved + nRemoved1;
    end
    return;
else
%     iClu1 = varargin{1};
%     P = varargin{2};
    nRemoved = 0;
    if ~isfield(P, 'nSkip_refrac'), P.nSkip_refrac = 4; end
    if ~isfield(P, 'fShow_refrac'), P.fShow_refrac = 1; end
    try
        viSpk1 = S_clu.cviSpk_clu{iClu1};
    catch
        viSpk1 = find(S_clu.viClu == iClu1);
    end
    if isempty(viSpk1), return; end

    viTime1 = viTime_spk(viSpk1);
    nRefrac = round(P.spkRefrac_ms * P.sRateHz / 1000);

    % removal loop
    vlKeep1 = true(size(viTime1));
    while (1)
        viKeep1 = find(vlKeep1);
        viRefrac11 = find(diff(viTime1(viKeep1)) < nRefrac) + 1;
        if isempty(viRefrac11), break; end

        vlKeep1(viKeep1(viRefrac11(1:P.nSkip_refrac:end))) = 0;
    end
    nRemoved = sum(~vlKeep1);
    nTotal1 = numel(vlKeep1);
    S_clu.viClu(viSpk1(~vlKeep1)) = 0;
    
    S_clu.cviSpk_clu{iClu1} = viSpk1(vlKeep1);
    S_clu.vnSpk_clu(iClu1) = sum(vlKeep1);
end

if get_(P, 'fVerbose')
    fprintf('Clu%d removed %d/%d (%0.1f%%) duplicate spikes\n', ...
        iClu1, nRemoved, nTotal1, nRemoved/nTotal1*100);
end
end %func


%--------------------------------------------------------------------------
function manual_test_(P, csCmd)
drawnow;
if nargin<2, csCmd = ''; end
if isempty(csCmd), csCmd = {'Mouse', 'Menu', 'FigWav', 'FigTime', 'FigWavCor', 'FigProj', 'Exit'}; end
if ischar(csCmd), csCmd = {csCmd}; end
S0 = get(0, 'UserData');
S_clu = S0.S_clu;
nClu = S0.S_clu.nClu;

for iCmd = 1:numel(csCmd)
    vcCmd1 = csCmd{iCmd};
    fprintf('\tTesting manual-mode %d/%d: %s\n', iCmd, numel(csCmd), vcCmd1);
    switch vcCmd1
        case 'Mouse' % simualte mouse click
            keyPress_fig_(get_fig_cache_('FigWav'), 'r');    %view whole
            fprintf('\tTesting mouse L/R clicks.\n');            
            viClu_test1 = [subsample_vr_(1:nClu, 5), nClu];
            for iClu1=viClu_test1
                fprintf('\t\tiCluCopy:%d/%d\n', iClu1, numel(viClu_test1));
                update_cursor_([], iClu1, 0); 
                keyPressFcn_cell_(get_fig_cache_('FigWav'), {'c','t','j','i','v','e','f'});
                drawnow;
                viClu_test2 = keep_lim_(iClu1 + [-2:2], [1, nClu]); 
                for iClu2=viClu_test2
                    fprintf('\t\t\tiCluPaste:%d/%d\n', iClu2, numel(viClu_test2));
                    update_cursor_([], iClu2, 1); 
                    keyPressFcn_cell_(get_fig_cache_('FigWav'), {'c','t','j','i','v','e','f'});
                    drawnow;
                end
            end
            
        case 'Menu' % run menu items, except for the exit and save (make a black list)
%             csMenu_skip = {'Show traces', 'Exit'};
%             hFigWav = get_fig_('FigWav');
            menu_test_(get_fig_('FigWav'), {'Show traces', 'Exit'});
%             vMenu0 = findobj('Type', 'uimenu', 'Parent', hFigWav);
%             cvMenu = cell(size(vMenu0));            
%             for iMenu0 = 1:numel(vMenu0)
%                 cvMenu{iMenu0} = findobj('Type', 'uimenu', 'Parent', vMenu0(iMenu0))';
%             end
%             vMenu = [cvMenu{:}];
%             cCallback_menu = get(vMenu, 'Callback');
%             csLabel_menu = get(vMenu, 'Label');
%             fprintf('\tTesting menu items\n');
%             for iMenu = 1:numel(csLabel_menu)    
%                 vcMenu = csLabel_menu{iMenu};             
%                 if ismember(vcMenu, csMenu_skip), continue; end
%                 try        
%                     hFunc = cCallback_menu{iMenu};
% %                     hFunc(hFigWav, []); %call function
%                     hFunc(vMenu(iMenu), []); %call function
%                     fprintf('\tMenu ''%s'' success.\n', vcMenu);                    
%                 catch
%                     fprintf(2, '\tMenu ''%s'' failed.\n', vcMenu);
%                     disperr_();
%                 end
%             end
            
        case 'FigWav' % test all possible keyboard press
            keyPress_fig_(get_fig_cache_('FigWav'), get_keyPress_('all'));          
            
        case 'FigTime'
            keyPress_fig_(get_fig_cache_('FigTime'), get_keyPress_('all'));          
            
        case 'FigWavCor'
            keyPress_fig_(get_fig_cache_('FigWavCor'), get_keyPress_('all'));           
            
        case 'FigProj'
            keyPress_fig_(get_fig_cache_('FigProj'), get_keyPress_('all'));                             
            
        case 'Exit'
            %fDebug_ui = 0;  set0_(fDebug_ui); % disable debug flag
            exit_manual_(get_fig_cache_('FigWav'));            
%             fDebug_ui = 1;  set0_(fDebug_ui);
            
        otherwise
            fprintf(2, 'Unsupported testing mode: %s\n', vcCmd1);
    end %swtich    
end
end %func


%--------------------------------------------------------------------------
function vrCentroid = com_(vrVpp, vrPos)
vrVpp_sq = vrVpp(:).^2;
vrCentroid = sum(vrVpp_sq .* vrPos(:)) ./ sum(vrVpp_sq);
end %func


%--------------------------------------------------------------------------
function viSpk_clu2 = spk_select_mid_(viSpk_clu1, P)
viTime_spk = get0_('viTime_spk');
iSpk_mid = round(numel(viTime_spk)/2);
viSpk_clu1_ord = rankorder_(abs(viSpk_clu1 - iSpk_mid), 'ascend');
nSpk1_max = round(numel(viSpk_clu1) / P.nTime_clu);
viSpk_clu2 = viSpk_clu1(viSpk_clu1_ord <= nSpk1_max);
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_update_(S_clu, iClu1, P)
% update cluster waveform and self correlation score
% mrWav not needed
global tnWav_spk tnWav_raw

% fSelect_mid = 1;

viSite_spk = get0_('viSite_spk');
% find clu center
viSpk_clu1 = find(S_clu.viClu == iClu1);
S_clu.cviSpk_clu{iClu1} = viSpk_clu1;
% if fSelect_mid  
%     iSite_clu1 = mode(viSite_spk(spk_select_mid_(viSpk_clu1, P)));
% else
% iSite_clu1 = mode(viSite_spk(viSpk_clu1));
% end
S_clu.viSite_clu(iClu1) = mode(viSite_spk(viSpk_clu1));
S_clu.vnSpk_clu(iClu1) = numel(viSpk_clu1);

% update mean waveform
S_clu = clu2wav_(S_clu, viSite_spk, tnWav_spk, tnWav_raw, iClu1);
vrSelfCorr_clu = get_diag_(S_clu.mrWavCor);
S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, iClu1);
S_clu.mrWavCor = set_diag_(S_clu.mrWavCor, vrSelfCorr_clu);
S_clu.mrWavCor(iClu1,iClu1) = S_clu_self_corr_(S_clu, tnWav_raw, iClu1);
S_clu = S_clu_position_(S_clu, iClu1);
S_clu = S_clu_quality_(S_clu, P, iClu1);
% [S_clu, S0] = S_clu_commit_(S_clu, 'S_clu_update_');
end %func


%--------------------------------------------------------------------------
function vhFig = get_fig_all_(csTag)
vhFig = nan(size(csTag));
for iFig=1:numel(csTag)
    try
        vhFig(iFig) = findobj('Tag', csTag{iFig}); 
    catch
%         disperr_();
    end
end %for
end %func


%--------------------------------------------------------------------------
function [hFig, S_fig] = get_fig_(vcTag, hFig)
% return figure handle based on the tag
% cache figure handles
% [usage]
% get_fig_(vcTag)
% get_fig_(vcTag, hFig) %build cache

% multiple tags requested
if iscell(vcTag), hFig = get_fig_all_(vcTag); return; end

S_fig = [];
try
    hFig = findobj('Tag', vcTag, 'Type', 'figure'); 
    if isempty(hFig)
        hFig = create_figure_(vcTag);
        try %set position if exists
            S0 = get(0, 'UserData');
            iFig = find(strcmp(S0.csFig, vcTag), 1, 'first');
            if ~isempty(iFig)
                set(hFig, 'OuterPosition', S0.cvrFigPos0{iFig});
            end
        catch
%             disperr_();
        end
    else
        hFig=hFig(end); %get later one
    end
    if nargout>1, S_fig = get(hFig, 'UserData'); end
catch
    hFig = [];
    disperr_();
end
end %end


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Type Figure added
function hFig = set_fig_(vcTag, S_fig)
% return figure handle based on the tag
% hFig = set_fig_(vcTag, S_fig)
% hFig = set_fig_(hFig, S_fig)
hFig = [];
try
    if ischar(vcTag)
        hFig = findobj('Tag', vcTag, 'Type', 'Figure');
    else
        hFig = vcTag;
    end
    set(hFig, 'UserData', S_fig); %figure property
catch
     disperr_();
end
end %end


%--------------------------------------------------------------------------
function keyPress_fig_(hFig, csKey)
% Simulate key press function
vcTag = get(hFig, 'Tag');
% S0 = get(0, 'UserData'); 
figure(hFig);
figure_wait_(1); 
event1.Key = '';
if ischar(csKey), csKey = {csKey}; end
nKeys = numel(csKey);
keyPressFcn__ = get(hFig, 'KeyPressFcn');
for i=1:nKeys
    try        
        event1.Key = csKey{i};
        keyPressFcn__(hFig, event1);
        fprintf('\tFigure ''%s'': Key ''%s'' success.\n', vcTag, csKey{i});
    catch
        fprintf(2, '\tFigure ''%s'': Key ''%s'' failed.\n', vcTag, csKey{i});
        disperr_();        
    end
%     pause(.1);
end
% drawnow;
figure_wait_(0);
end %func


%--------------------------------------------------------------------------
function csKeys = get_keyPress_(vcType)
% return key press
switch vcType
    case 'all'
        csKeys = [get_keyPress_('arrows'), get_keyPress_('misc'), get_keyPress_('alphanumeric')];
    case 'alphanumeric'
        csKeys = char([double('0'):double('9'), double('A'):double('Z'), double('a'):double('z')]);
        csKeys = num2cell(csKeys);
    case 'arrows'
        csKeys = {'uparrow', 'downarrow', 'leftarrow', 'rightarrow'};
    case 'misc'
        csKeys = {'home', 'end', 'space', 'esc'};
end %switch
end %func


%--------------------------------------------------------------------------
function hPoly = impoly_(varargin)
global fDebug_ui
% if get_set_([], 'fDebug_ui', 0)
if fDebug_ui==1
    hPoly = []; %skip the test if debugging
else
    hPoly = impoly(varargin{:});
end
end


%--------------------------------------------------------------------------
function vi = keep_lim_(vi, lim)
vi = vi(vi>=lim(1) & vi <= lim(end));
end


%--------------------------------------------------------------------------
function [S_fig, maxAmp_prev, hFig] = set_fig_maxAmp_(vcFig, event)
[hFig, S_fig] = get_fig_cache_(vcFig);
if isempty(S_fig)
    P = get0_('P');
    S_fig.maxAmp = P.maxAmp;
end
maxAmp_prev = S_fig.maxAmp;
S_fig.maxAmp = change_amp_(event, maxAmp_prev);
set(hFig, 'UserData', S_fig);
end


%--------------------------------------------------------------------------
function S = rmfield_(S, varargin)
% varargin: list of fields to remove
for i=1:numel(varargin)
    if isfield(S, varargin{i})
        S = rmfield(S, varargin{i});
    end
end
end %func


%--------------------------------------------------------------------------
function hFig = create_figure_(vcTag, vrPos, vcName, fToolbar, fMenubar)
% or call external create_figure()
if nargin<2, vrPos = []; end
if nargin<3, vcName = ''; end
if nargin<4, fToolbar = 0; end
if nargin<5, fMenubar = 0; end
if isempty(vcTag)
    hFig = figure();
else
    hFig = figure_new_(vcTag); 
end
set(hFig, 'Name', vcName, 'NumberTitle', 'off', 'Color', 'w');
clf(hFig);
set(hFig, 'UserData', []); %empty out the user data
if ~fToolbar
    set(hFig, 'ToolBar', 'none'); 
else
    set(hFig, 'ToolBar', 'figure'); 
end
if ~fMenubar
    set(hFig, 'MenuBar', 'none'); 
else
    set(hFig, 'MenuBar', 'figure'); 
end

if ~isempty(vrPos), resize_figure_(hFig, vrPos); end
clf(hFig);
end %func


%--------------------------------------------------------------------------
function hAx = axes_new_(hFig)
if ischar(hFig), hFig = get_fig_(hFig); end
figure(hFig); %set focus to figure %might be slow
clf(hFig); 
hAx = axes(); 
hold(hAx, 'on');
end %func


%--------------------------------------------------------------------------
function hFig = figure_new_(vcTag)
%remove prev tag duplication
delete_multi_(findobj('Tag', vcTag, 'Type', 'Figure')); 

hFig = figure('Tag', vcTag);
end %func


%--------------------------------------------------------------------------
function hFig = resize_figure_(hFig, posvec0, fRefocus)
if nargin<3, fRefocus = 1; end
height_taskbar = 40;

pos0 = get(groot, 'ScreenSize'); 
width = pos0(3); 
height = pos0(4) - height_taskbar;
% width = width;
% height = height - 132; %width offset
% width = width - 32;
posvec = [0 0 0 0];
posvec(1) = max(round(posvec0(1)*width),1);
posvec(2) = max(round(posvec0(2)*height),1) + height_taskbar;
posvec(3) = min(round(posvec0(3)*width), width);
posvec(4) = min(round(posvec0(4)*height), height);
% drawnow;
if isempty(hFig)
    hFig = figure; %create a figure
else
    hFig = figure(hFig);
end
drawnow;
set(hFig, 'OuterPosition', posvec, 'Color', 'w', 'NumberTitle', 'off');
end


%--------------------------------------------------------------------------
function traces_test_(P)
drawnow;
% csCmd = {'Mouse', 'Menu', 'FigWav', 'FigTime', 'FigWavCor', 'FigProj', 'Exit'};

% for iCmd = 1:numel(csCmd)
% vcCmd1 = csCmd{iCmd};
% fprintf('\tTesting manual-mode %d/%d: %s\n', iCmd, numel(csCmd), vcCmd1);
hFig = get_fig_('Fig_traces');
keyPress_fig_(hFig, get_keyPress_('all'));
try
    close(hFig); %close traces figure. other figures may remain
    close(get_fig_('FigPsd'));
catch
end
end %func


%--------------------------------------------------------------------------
function gui_test_(P, vcFig, csMenu_skip)
if nargin<3, csMenu_skip = {}; end
drawnow;
hFig = get_fig_(vcFig);
keyPress_fig_(hFig, get_keyPress_('all'));

% Menu test
menu_test_(hFig, csMenu_skip);

try
    close(hFig); %close traces figure. other figures may remain
catch
    ;
end
end %func


%--------------------------------------------------------------------------
function [vlSuccess_menu, csLabel_menu] = menu_test_(hFig, csMenu_skip)
vMenu0 = findobj('Type', 'uimenu', 'Parent', hFig);
cvMenu = cell(size(vMenu0));            
for iMenu0 = 1:numel(vMenu0)
    cvMenu{iMenu0} = findobj('Type', 'uimenu', 'Parent', vMenu0(iMenu0))';
end
vMenu = [cvMenu{:}];
cCallback_menu = get(vMenu, 'Callback');
csLabel_menu = get(vMenu, 'Label');
fprintf('\tTesting menu items\n');

vlSuccess_menu = true(size(csLabel_menu));
for iMenu = 1:numel(csLabel_menu)    
    vcMenu = csLabel_menu{iMenu};             
    if ismember(vcMenu, csMenu_skip), continue; end
    try        
        hFunc = cCallback_menu{iMenu};
        if isempty(hFunc), continue; end
%                     hFunc(hFigWav, []); %call function
        hFunc(vMenu(iMenu), []); %call function
        fprintf('\tMenu ''%s'' success.\n', vcMenu);             
    catch
        fprintf(2, '\tMenu ''%s'' failed.\n', vcMenu);
        disperr_();
        vlSuccess_menu(iMenu) = 0;
    end
end
end %func


%--------------------------------------------------------------------------
function hRect = imrect_(varargin)
global fDebug_ui
    
hRect = []; %skip the test if debugging
% if get_set_([], 'fDebug_ui', 0) && nargin < 2
if fDebug_ui==1 && nargin < 2
    return;
else
    try
        hRect = imrect(varargin{:});
    catch
        fprintf(2, 'Install image processing toolbox\n');
    end
end
end


%--------------------------------------------------------------------------
function csAns = inputdlg_(varargin)
% return default answer
global fDebug_ui
% if get_set_([], 'fDebug_ui', 0)
if fDebug_ui==1
    if numel(varargin)==4
        csAns = varargin{4};
    else
        csAns = [];
    end
else
    csAns = inputdlg(varargin{:});
end
end


%--------------------------------------------------------------------------
function vr_uV = bit2uV_(vn, P)
% use only for filtered traces

if nargin<2, P = get0_('P'); end
% if isempty(get_(P, 'nDiff_filt')), P.nDiff_filt = 0; end
if strcmpi(get_filter_(P), 'sgdiff')
    norm = sum((1:P.nDiff_filt).^2) * 2;
else
    norm = 1;
end
% switch P.nDiff_filt
%     case 0, norm = 1;
%     case 1, norm = 2;
%     case 2, norm = 10;
%     case 3, norm = 28;
%     otherwise, norm = 60;
% end %switch
vr_uV = single(vn) * single(P.uV_per_bit / norm);
end


%--------------------------------------------------------------------------
function vr_uV = uV2bit_(vn, P)
% use only for filtered traces

if nargin<2, P = get0_('P'); end
if isempty(P.nDiff_filt), P.nDiff_filt = 0; end
switch P.nDiff_filt
    case 0, norm = 1;
    case 1, norm = 2;
    case 2, norm = 10;
    case 3, norm = 28;
    otherwise, norm = 60;
end %switch
vr_uV = single(vn) / single(P.uV_per_bit / norm);
end


%--------------------------------------------------------------------------
function flag = isvalid_(h)
try
    flag = isvalid(h);
catch
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Generalized to any figures previously querried.
function [hFig, S_fig] = get_fig_cache_(vcFig_tag)
% clear before starting manual
% Return from persistent cache
% Create a new figure if Tag doesn't exist

persistent S_fig_cache_
% persistent hFigPos hFigMap hFigWav hFigTime hFigProj hFigWavCor hFigHist hFigIsi hFigCorr hFigRD hFig_traces hFig_preview
if iscell(vcFig_tag)
    if nargout==1
        hFig = cellfun(@(vc)get_fig_cache_(vc), vcFig_tag, 'UniformOutput', 0);
    else
        [hFig, S_fig] = cellfun(@(vc)get_fig_cache_(vc), vcFig_tag, 'UniformOutput', 0);
    end
    return; 
end
if isempty(S_fig_cache_)
    hFig = get_fig_(vcFig_tag);
    S_fig_cache_ = struct(vcFig_tag, hFig);
else
    if isfield(S_fig_cache_, vcFig_tag)
        hFig = S_fig_cache_.(vcFig_tag);
        if isvalid_(hFig)
            hFig = S_fig_cache_.(vcFig_tag);            
        else
            hFig = get_fig_(vcFig_tag);
            S_fig_cache_.(vcFig_tag) = hFig;
        end
    else
        hFig = get_fig_(vcFig_tag);
        S_fig_cache_.(vcFig_tag) = hFig;
    end
end
% switch vcFig
%     case 'FigPos'
%         if ~isvalid_(hFigPos), hFigPos = get_fig_(vcFig); end
%         hFig = hFigPos;    
%     case 'FigMap'
%         if ~isvalid_(hFigMap), hFigMap = get_fig_(vcFig); end
%         hFig = hFigMap;        
%     case 'FigWav'
%         if ~isvalid_(hFigWav), hFigWav = get_fig_(vcFig); end
%         hFig = hFigWav;     
%     case 'FigTime'
%         if ~isvalid_(hFigTime), hFigTime = get_fig_(vcFig); end
%         hFig = hFigTime;                
%     case 'FigProj'
%         if ~isvalid_(hFigProj), hFigProj = get_fig_(vcFig); end
%         hFig = hFigProj;
%     case 'FigWavCor'
%         if ~isvalid_(hFigWavCor), hFigWavCor = get_fig_(vcFig); end
%         hFig = hFigWavCor;        
%     case 'FigHist'
%         if ~isvalid_(hFigHist), hFigHist = get_fig_(vcFig); end
%         hFig = hFigHist;                
%     case 'FigIsi'
%         if ~isvalid_(hFigIsi), hFigIsi = get_fig_(vcFig); end
%         hFig = hFigIsi;   
%     case 'FigCorr'
%         if ~isvalid_(hFigCorr), hFigCorr = get_fig_(vcFig); end
%         hFig = hFigCorr;           
%     case 'FigRD'
%         if ~isvalid_(hFigRD), hFigRD = get_fig_(vcFig); end
%         hFig = hFigRD;   
%     case 'Fig_traces'
%         if ~isvalid_(hFig_traces), hFig_traces = get_fig_(vcFig); end
%         hFig = hFig_traces;      
%     case 'Fig_preview'
%         if ~isvalid_(hFig_preview), hFig_preview = get_fig_(vcFig); end
%         hFig = hFig_preview;            
%     otherwise
%         fprintf(2, 'get_fig_cache_: invalid figure tag: %s\n', vcFig); 
% end %switch
if nargout>1, S_fig = get(hFig, 'UserData'); end
end %func


%--------------------------------------------------------------------------
function import_jrc1_(vcFile_prm)
% import jrc1 so that i can visualize the output
global tnWav_raw tnWav_spk trFet_spk
% convert jrc1 format (_clu and _evt) to jrc3 format. no overwriting
% receive spike location, time and cluster number. the rest should be taken care by jrc3 processing

% Load info from previous version: time, site, spike
P = loadParam_(vcFile_prm);
Sevt = load(strrep(P.vcFile_prm, '.prm', '_evt.mat'));
if isfield(Sevt, 'Sevt'), Sevt = Sevt.Sevt; end
S0 = struct('viTime_spk', Sevt.viSpk, 'viSite_spk', Sevt.viSite, 'P', P);

[tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, S0.viTime_spk, S0.viSite_spk);
set(0, 'UserData', S0);

% Save to file
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk);
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk);
save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));

% cluster and describe
sort_(P);
S0 = get(0, 'UserData');
Sclu = load_(strrep(P.vcFile_prm, '.prm', '_clu.mat'));
if isfield(Sclu, 'Sclu'), Sclu = Sclu.Sclu; end
if ~isempty(Sclu)
    S0.S_clu.viClu = Sclu.viClu; %skip FigRD step for imported cluster
end
set(0, 'UserData', S0);

describe_(P.vcFile_prm);
end %func


%--------------------------------------------------------------------------
function S_mat = load_(vcFile, csVar, fVerbose)
% return empty if the file doesn't exist
% load_(vcFile)
% load_(vcFile, vcVar)
% load_(vcFile, csVar)
% load_(vcFile, csVar, fVerbose)
% load_(vcFile, [], fVerbose)

S_mat = [];
if nargin<2, csVar={}; end
if ischar(csVar), csVar = {csVar}; end
if nargin<3, fVerbose=1; end
if exist(vcFile, 'file')
    try
        S_mat = load(vcFile);
    catch
        fprintf(2, 'Invalid .mat format: %s\n', vcFile);
    end
else
    if fVerbose
        fprintf(2, 'File does not exist: %s\n', vcFile);
    end
end
if ~isempty(csVar)
    S_mat = get_(S_mat, csVar{:});
end
end %func


%--------------------------------------------------------------------------
function [viTime_spk, vrAmp_spk, viSite_spk] = detect_spikes_(mnWav3, vnThresh_site, vlKeep_ref, P)
% fMerge_spk = 1;
fMerge_spk = get_set_(P, 'fMerge_spk', 1);

[n1, nSites, ~] = size(mnWav3);
[cviSpk_site, cvrSpk_site] = deal(cell(nSites,1));

fprintf('\tDetecting spikes from each channel.\n\t\t'); t1=tic;
% parfor iSite = 1:nSites   
for iSite = 1:nSites   
    % Find spikes
    [viSpk11, vrSpk11] = spikeDetectSingle_fast_(mnWav3(:,iSite), P, vnThresh_site(iSite));
    fprintf('.');
    
    % Reject global mean
    if isempty(vlKeep_ref)
        cviSpk_site{iSite} = viSpk11;
        cvrSpk_site{iSite} = vrSpk11;        
    else
        [cviSpk_site{iSite}, cvrSpk_site{iSite}] = select_vr_(viSpk11, vrSpk11, find(vlKeep_ref(viSpk11)));
    end
end
vnThresh_site = gather_(vnThresh_site);
nSpks1 = sum(cellfun(@numel, cviSpk_site));
fprintf('\n\t\tDetected %d spikes from %d sites; took %0.1fs.\n', nSpks1, nSites, toc(t1));

% Group spiking events using vrWav_mean1. already sorted by time
if fMerge_spk
    fprintf('\tMerging spikes...'); t2=tic;
    [viTime_spk, vrAmp_spk, viSite_spk] = spikeMerge_(cviSpk_site, cvrSpk_site, P);
    fprintf('\t%d spiking events found; took %0.1fs\n', numel(viSite_spk), toc(t2));
else
    viTime_spk = cell2mat_(cviSpk_site);
    vrAmp_spk = cell2mat_(cvrSpk_site);
    viSite_spk = cell2vi_(cviSpk_site);
    %sort by time
    [viTime_spk, viSrt] = sort(viTime_spk, 'ascend');
    [vrAmp_spk, viSite_spk] = multifun_(@(x)x(viSrt), vrAmp_spk, viSite_spk);
end
vrAmp_spk = gather_(vrAmp_spk);
end %func


%--------------------------------------------------------------------------
function [viTime_spk11, viSite_spk11] = filter_spikes_(viTime_spk0, viSite_spk0, tlim)
% Filter spikes that is within tlim specified

[viTime_spk11, viSite_spk11] = deal([]); 
if isempty(viTime_spk0), return; end
viKeep11 = find(viTime_spk0 >= tlim(1) & viTime_spk0 <= tlim(end));
viTime_spk11 = viTime_spk0(viKeep11)  + (1 - tlim(1)); % shift spike timing
if ~isempty(viSite_spk0)
    viSite_spk11 = viSite_spk0(viKeep11);
else
    viSite_spk11 = [];
end
end %func


%--------------------------------------------------------------------------
function import_silico_(vcFile_prm, fSort)
% need _gt struct?
% import silico ground truth
% S_gt: must contain viTime, viSite, viClu 
% [usage]
% import_silico_(vcFile_prm, 0): use imported sort result (default)
% import_silico_(vcFile_prm, 1): use jrclust sort result
if nargin<2, fSort = 0; end

global tnWav_raw tnWav_spk trFet_spk
% convert jrc1 format (_clu and _evt) to jrc3 format. no overwriting
% receive spike location, time and cluster number. the rest should be taken care by jrc3 processing
P = loadParam_(vcFile_prm); %makeParam_kilosort_
if isempty(P), return; end

S_gt = load(strrep(vcFile_prm, '.prm', '_gt.mat'));   %must contain viTime, viSite, viClu
% vnSpk = cellfun(@numel, S.a);    
% viClu = int32(cell2mat_(arrayfun(@(n)n*ones(1, vnSpk(n)), 1:numel(vnSpk), 'UniformOutput', 0)));
% viTime = int32(cell2mat_(S.a) * 20); % Convert to sample # (saved in ms unit & sampling rate =20KHZ)
if ~isfield(S_gt, 'viSite')
    S_gt.viSite = S_gt.viSite_clu(S_gt.viClu); 
end
S0 = struct('viTime_spk', S_gt.viTime(:), 'viSite_spk', S_gt.viSite(:), 'P', P, 'S_gt', S_gt);

[tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, S0.viTime_spk, S0.viSite_spk);
set(0, 'UserData', S0);

% Save to file
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk); 
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk);

% cluster and describe
S0 = sort_(P);
if ~fSort %use ground truth cluster 
    S0.S_clu = S_clu_new_(S_gt.viClu, S0);    
end
set(0, 'UserData', S0);

% Save
save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
describe_(S0);
end %func


%--------------------------------------------------------------------------
function import_ksort_(vcFile_prm, fSort)
% import_ksort_(P, fSort)
% import_ksort_(vcFile_prm, fSort)
fMerge_post = 0;
% import kilosort result
% fSort: do sort using jrc3
if isstruct(vcFile_prm)
    P = vcFile_prm;
    vcFile_prm = P.vcFile_prm;
else
    P = loadParam_(vcFile_prm); %makeParam_kilosort_
end

% S_ksort = load(strrep(P.vcFile_prm, '.prm', '_ksort.mat')); % contains rez structure
if nargin<2, fSort = 0; end %
global tnWav_raw tnWav_spk trFet_spk
% convert jrc1 format (_clu and _evt) to jrc3 format. no overwriting
% receive spike location, time and cluster number. the rest should be taken care by jrc3 processing

% Create a prm file to start with. set the filter parameter correctly. features? 
if isempty(P), return; end
S_ksort = load(strrep(P.vcFile_prm, '.prm', '_ksort.mat')); %get site # and 
viTime_spk = S_ksort.rez.st3(:,1) - 6; %spike time (apply shift factor)
viClu = S_ksort.rez.st3(:,2); % cluster

viClu_post = 1 + S_ksort.rez.st3(:,5);  %post-merging result
tnWav_clu = S_ksort.rez.Wraw; %nC, nT, nClu
tnWav_clu = -abs(tnWav_clu);
tnWav_clu = permute(tnWav_clu, [2,1,3]);
mnMin_clu = squeeze_(min(tnWav_clu, [], 1));
[~, viSite_clu] = min(mnMin_clu, [], 1); %cluster location
viSite = 1:numel(P.viSite2Chan);
viSite(P.viSiteZero) = [];
viSite_clu = viSite(viSite_clu);
viSite_spk = viSite_clu(viClu);
% vnAmp_spk = S_ksort.rez.st3(:,3);

S0 = struct('viTime_spk', int32(viTime_spk), 'viSite_spk', int32(viSite_spk), 'P', P, 'S_ksort', S_ksort);
[tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, S0.viTime_spk, S0.viSite_spk);
set(0, 'UserData', S0);

% Save to file
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkraw.jrc'), tnWav_raw);
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkwav.jrc'), tnWav_spk); 
write_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), trFet_spk); 

% cluster and describe
S0 = sort_(P);
if ~fSort %use ground truth cluster 
    if fMerge_post
        S0.S_clu = S_clu_new_(viClu_post, S0);    
    else
        S0.S_clu = S_clu_new_(viClu, S0);    
    end
    S0.S_clu = S_clu_sort_(S0.S_clu, 'viSite_clu'); 
end
S0.S_clu = S_clu_new_(S0.S_clu);
set(0, 'UserData', S0);

% Save
save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
describe_(S0);
end %func


%--------------------------------------------------------------------------
function vl = matchFileExt_(csFiles, vcExt, vlDir)
% vcExt can be a cell
% ignore dir
% matchFileExt_(csFiles, vcExt, vlDir)
% matchFileExt_(csFiles, csExt, vlDir) %multiple extension check
if ischar(csFiles), csFiles={csFiles}; end
vl = false(size(csFiles));

for i=1:numel(csFiles)
    [~,~,vcExt1] = fileparts(csFiles{i});
    vl(i) = any(strcmpi(vcExt1, vcExt));
end
if nargin >= 3
    vl = vl | vlDir; %matches if it's directory
end
end %func


%--------------------------------------------------------------------------
function P = file2struct_(vcFile_file2struct)
% Run a text file as .m script and result saved to a struct P
% _prm and _prb can now be called .prm and .prb files

try 
    P = file2struct__(vcFile_file2struct); % new version
catch
    P = file2struct_1_(vcFile_file2struct); % old version
end
% if isempty(P)
%     copyfile(vcFile_file2struct, 'temp_eval.m', 'f');
%     try
%         eval('temp_eval.m');
%     catch
%         disp(lasterr());
%     end
% end
end %func


%--------------------------------------------------------------------------
function P = file2struct_1_(vcFile_file2struct)
if ~exist(vcFile_file2struct, 'file')
    fprintf(2, '%s does not exist.\n', vcFile_file2struct);
    P = [];
    return;
end

% load text file
fid=fopen(vcFile_file2struct, 'r');
csCmd = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
csCmd = csCmd{1};

% parse command
for iCmd=1:numel(csCmd)
    try
        vcLine1 = strtrim(csCmd{iCmd});
        if isempty(vcLine1), continue; end
        if find(vcLine1=='%', 1, 'first')==1, continue; end
        iA = find(vcLine1=='=', 1, 'first');
        if isempty(iA), continue; end            
        iB = find(vcLine1=='(', 1, 'first');
        if ~isempty(iB) && iB<iA, iA=iB; end
        eval(vcLine1);
        vcVar1 = strtrim(vcLine1(1:iA-1));
        eval(sprintf('P.(vcVar1) = %s;', vcVar1));
    catch
        fprintf(2, lasterr);
    end
end %for
end %func


%--------------------------------------------------------------------------
% function P = appendStruct_(P, varargin)
% % backward compatibility
% P = struct_merge_(P, varargin{:});
% end


%--------------------------------------------------------------------------
function vcPath = replacePath_(vcPath1, vcPath2)
% replace path1 with path2
[~, vcFname1, vcExt1] = fileparts(vcPath1);
[vcDir2,~,~] = fileparts(vcPath2);
if ~isempty(vcDir2)
    vcPath = [vcDir2, filesep(), vcFname1, vcExt1];
else
    vcPath = [vcFname1, vcExt1];
end
end %func


%---------------------------------------------------------------------------
function out = ifeq_(if_, true_, false_)
if (if_)
    out = true_;
else
    out = false_;
end
end %func


%--------------------------------------------------------------------------
function figure_wait_(fWait, vhFig)
% set all figures pointers to watch
if nargin<2
    vhFig = gcf;
%     csFig = get0_('csFig');
%     if isempty(csFig)
%         vhFig = findobj('Type', 'Figure');
%     else
%         vhFig = get_fig_cache_(csFig);
%     end
end
if fWait    
    set_(vhFig, 'Pointer', 'watch');
else
    set_(vhFig, 'Pointer', 'arrow');
end
end %func


%--------------------------------------------------------------------------
function nBytes = getBytes_(vcFile)
S_dir = dir(vcFile);
if isempty(S_dir), nBytes=[]; return; end
nBytes = S_dir(1).bytes;
end %func


%--------------------------------------------------------------------------
function S = makeStruct_(varargin)
%MAKESTRUCT all the inputs must be a variable. 
%don't pass function of variables. ie: abs(X)
%instead create a var AbsX an dpass that name

S=[];
for i=1:nargin
    S = setfield(S, inputname(i), varargin{i});
end
end %func


%--------------------------------------------------------------------------
function [mr, vi] = subsample_mr_(mr, nMax, dimm)
%[mr, vi] = subsample_mr_(mr, nMax, dimm)
% subsample the column
if nargin<3, dimm = 2; end
if isempty(nMax), return ;end

n = size(mr,dimm);
nSkip = max(floor(n / nMax), 1);
vi = 1:nSkip:n;
if nSkip==1, return; end
vi = vi(1:nMax);

switch dimm
    case 2
        mr = mr(:,vi);
    case 1
        mr = mr(vi,:);
end

if nargout>=2
    if n > nMax
        vi = 1:nSkip:n;
        vi = vi(1:nMax);
    else
        vi = 1:n;
    end
end
end %func


%--------------------------------------------------------------------------
function S_old = struct_append_(S_old, S_new, varargin)
if isempty(S_new), return; end
if nargin==2
    varargin = fieldnames(S_new);
end   
for i=1:numel(varargin)
    try
        S_old.(varargin{i}) = S_new.(varargin{i});
    catch
        ;
    end
end %for
end %func


%--------------------------------------------------------------------------
function varargout = multifun_(hFun, varargin)
% apply same function to the input, unary function only

if nargout ~= numel(varargin), error('n arg mismatch'); end
for i=1:nargout
    try
        varargout{i} = hFun(varargin{i});
    catch
        varargout{i} = varargin{i};
    end
end
end %func


%--------------------------------------------------------------------------
function mr = reshape_vr2mr_(vr, nwin)
nbins = ceil(numel(vr)/nwin);
vr(nbins*nwin) = 0; %expand size
mr = reshape(vr(1:nbins*nwin), nwin, nbins);
end %func


%--------------------------------------------------------------------------
% 9/29/17 Updating the version number when saving JRCLUST
function S0 = save0_(vcFile_mat)
% save S0 structure to a mat file
try
    fprintf('Saving 0.UserData to %s...\n', vcFile_mat);
    warning off;
    S0 = get(0, 'UserData'); %add gather script
    
    % update version number
    S0.P.version = jrc_version_();
    P = S0.P;
    set0_(P);
    
    struct_save_(S0, vcFile_mat, 1);
    vcFile_prm = S0.P.vcFile_prm;
    export_prm_(vcFile_prm, strrep(vcFile_prm, '.prm', '_full.prm'), 0);
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function [vl, vr] = thresh_mad_(vr, thresh_mad)
% single sided, no absolute value

nsubs = 300000;
offset = median(subsample_vr_(vr, nsubs));
vr = vr - offset; %center the mean
factor = median(abs(subsample_vr_(vr, nsubs)));
if isempty(thresh_mad) || thresh_mad==0
    vl = true(size(vr));
else
    vl = vr < factor * thresh_mad;
end
if nargout>=2
    vr = vr / factor; %MAD unit
end
end %func


%--------------------------------------------------------------------------
function vi = subsample_vr_(vi, nMax)
if numel(vi)>nMax
    nSkip = floor(numel(vi)/nMax);
    if nSkip>1, vi = vi(1:nSkip:end); end
    if numel(vi)>nMax, vi = vi(1:nMax); end
end
end %func


%--------------------------------------------------------------------------
function [vrPow, vrFreq] = plotMedPower_(mrData, varargin)

if numel(varargin) ==1
    if ~isstruct(varargin{1}), P.sRateHz = varargin{1}; 
    else P = varargin{1};
    end
else
    P = funcInStr_(varargin{:});
end
    
P = funcDefStr_(P, 'viChanExcl', [1 18 33 50 65 82 97 114], 'sRateHz', 25000, ...
    'nSmooth', 3, 'LineStyle', 'k-', 'fPlot', 1, 'fKHz', 0, 'vcMode', 'max');
    
if iscell(mrData) %batch mode
    csFname = mrData;
    [vrPow, vrFreq] =deal(cell(numel(csFname)));
    for i=1:numel(csFname)
        hold on;
        [vrPow{i}, vrFreq{i}] = plotMedPower_(csFname{i});
    end
    legend(csFname);
    return;
else
    vcFname='';
end
warning off;

mrData = fft(mrData);
mrData = real(mrData .* conj(mrData)) / size(mrData,1) / (P.sRateHz/2);
% mrPowFilt = filter([1 1 1], 3, mrPow);
% mrPow = fftshift(mrPow);
imid0 = ceil(size(mrData,1)/2);
vrFreq = (0:size(mrData,1)-1) * (P.sRateHz/size(mrData,1));
vrFreq = vrFreq(2:imid0);
viChan = setdiff(1:size(mrData,2), P.viChanExcl);
if size(mrData,2)>1
    switch P.vcMode
        case 'mean'
            vrPow = mean(mrData(2:imid0,viChan), 2);    
        case 'max'
            vrPow = max(mrData(2:imid0,viChan), [], 2);    
    end
else
    vrPow = mrData(2:imid0,1);
end
% vrPow = std(mrData(:,viChan), 1, 2);
if P.nSmooth>1, vrPow = filterq_(ones([P.nSmooth,1]),P.nSmooth,vrPow); end

if P.fPlot
    if P.fKHz, vrFreq = vrFreq/1000; end
    plot(vrFreq, pow2db_(vrPow), P.LineStyle);
%     set(gca, 'YScale', 'log'); 
%     set(gca, 'XScale', 'linear'); 
    xlabel('Frequency (Hz)'); ylabel('Mean power across sites (dB uV^2/Hz)');
    % xlim([0 P.sRateHz/2]);
    grid on;
    try
    xlim(vrFreq([1, end]));
    set(gcf,'color','w');
    title(vcFname, 'Interpreter', 'none');
    catch
        
    end
end
end %func


%--------------------------------------------------------------------------
function P = funcInStr_( varargin )

if isempty(varargin), P=struct(); return; end
if isstruct(varargin{1}), P = varargin{1}; return; end

csNames = varargin(1:2:end);
csValues = varargin(2:2:end);
P = struct();
for iField=1:numel(csNames)
    if ~isfield(P, csNames{iField})
%         v1 = csValues{iField};
%         eval(sprintf('P.%s = v1;', csNames{iField}));
        P = setfield(P, csNames{iField}, csValues{iField});
    end
end
end %func


%--------------------------------------------------------------------------
function P = funcDefStr_(P, varargin)
csNames = varargin(1:2:end);
csValues = varargin(2:2:end);

for iField=1:numel(csNames)
    if ~isfield(P, csNames{iField})
        P = setfield(P, csNames{iField}, csValues{iField});
    end
end
end %func


%--------------------------------------------------------------------------
function varargout = select_vr_(varargin)
% [var1, var2, ...] = select_vr(var1, var2, ..., index)

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


%--------------------------------------------------------------------------
function S0 = load0_(vcFile_mat)
% Load a mat file structure and set to 0 structure
% S0 = load0_(vcFile_mat)
% S0 = load0_(P)
% only set the S0 if it's found
if isstruct(vcFile_mat)
    P = vcFile_mat;
    vcFile_mat = strrep(P.vcFile_prm, '.prm', '_jrc.mat');
end

if ~exist(vcFile_mat, 'file')
    S0 = [];
    fprintf('File %s does not exist\n', vcFile_mat);
    return;
end
    
try
    fprintf('loading %s...\n', vcFile_mat); t1=tic;
    S0 = load(vcFile_mat);    
    set(0, 'UserData', S0);
    fprintf('\ttook %0.1fs\n', toc(t1));
catch
    S0 = [];
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function S = struct_add_(S, varargin)

for i=1:numel(varargin)
    S.(inputname(i+1)) = varargin{i};
end
end %func


%--------------------------------------------------------------------------
function [tr, miRange] = mr2tr3_(mr, spkLim, viTime, viSite, fMeanSubt)
% tr: nSamples x nSpikes x nChans

if nargin<4, viSite=[]; end %faster indexing
if nargin<5, fMeanSubt=0; end

% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype
if isempty(viTime), tr=[]; return; end
[N, M] = size(mr);
if ~isempty(viSite), M = numel(viSite); end
if iscolumn(viTime), viTime = viTime'; end

viTime0 = [spkLim(1):spkLim(end)]'; %column
miRange = bsxfun(@plus, int32(viTime0), int32(viTime));
miRange = min(max(miRange, 1), N);
miRange = miRange(:);

if isempty(viSite)
    tr = mr(miRange,:);
else
    tr = mr(miRange, viSite);
end
tr = reshape(tr, [numel(viTime0), numel(viTime), M]);

if fMeanSubt
%     trWav1 = single(permute(trWav1, [1,3,2])); 
    tr = single(tr);
    dimm1 = size(tr);
    tr = reshape(tr, size(tr,1), []);
    tr = bsxfun(@minus, tr, mean(tr)); %mean subtract
    tr = reshape(tr, dimm1);    
end
end %func


%--------------------------------------------------------------------------
function disp_stats_(vr)
% fprintf('n, mu/sd, (med,q25,q75), min-max:\t %d, %0.2f/%0.2f, (%0.2f, %0.2f, %0.2f), %0.2f-%0.2f\n', ...
%     numel(vr), mean(vr), std(vr), quantile(vr, [.5, .25, .75]), min(vr), max(vr));
fprintf('n, mu/sd, (10,25,*50,75,90%%), min-max:\t %d, %0.2f/%0.2f, (%0.2f, %0.2f, *%0.2f, %0.2f, %0.2f), %0.2f-%0.2f\n', ...
    numel(vr), mean(vr), std(vr), quantile(vr, [.1,.25,.5,.75,.9]), min(vr), max(vr));
end %func

%--------------------------------------------------------------------------
function [mrWav, S_tsf] = importTSF_(fname, varargin)
fid = fopen(fname, 'r');
S_tsf = importTSF_header_(fid);
n_vd_samples = S_tsf.n_vd_samples;
n_electrodes = S_tsf.nChans;
mrWav = reshape(fread(fid, n_vd_samples * n_electrodes, '*int16'), [n_vd_samples,n_electrodes]);
fclose(fid);
end %func


%--------------------------------------------------------------------------
function Sfile = importTSF_header_(arg1)
% Sfile = importTSF_header(vcFname) %pass file name string
% Sfile = importTSF_header(fid) %pass fid

if ischar(arg1)
    fid = fopen(arg1, 'r');
    vcFname = arg1;
else
    fid = arg1;
    vcFname = [];
end

% 16 bytes
header = fread(fid, 16, 'char*1=>char');

% 4x5 bytes
iformat = fread(fid, 1, 'int32');
SampleFrequency = fread(fid, 1, 'int32');
n_electrodes = fread(fid, 1, 'int32');
n_vd_samples = fread(fid, 1, 'int32');
vscale_HP = fread(fid, 1, 'single');

% electrode info
Siteloc = zeros(2, n_electrodes, 'int16'); %2 x n_elec
Readloc = zeros(1, n_electrodes, 'int32'); %read location
for i_electrode = 1:n_electrodes
    Siteloc(:, i_electrode) = fread(fid, 2, 'int16');
    Readloc(i_electrode) = fread(fid, 1, 'int32');
end
if ~isempty(vcFname), fclose(fid); end

Sfile = struct('header', header, 'iformat', iformat, ...
    'sRateHz', SampleFrequency, 'nChans', n_electrodes, ...
    'n_vd_samples', n_vd_samples, 'vscale_HP', vscale_HP, ...
    'Siteloc', Siteloc, 'tLoaded', n_vd_samples/SampleFrequency, 'Readloc', Readloc);
end %func


%--------------------------------------------------------------------------
function dc = calc_dc2_(S0, P, vlRedo_spk)
global trFet_spk
if nargin<3, vlRedo_spk=[]; end
fprintf('Calculating Dc\n\t'); t1=tic;
nSites = numel(P.viSite2Chan);
vrDc2_site = nan(1, nSites);
for iSite = 1:nSites
    [mrFet12_, viSpk12_, n1_, n2_, viiSpk12_ord_] = fet12_site_(trFet_spk, S0, P, iSite, vlRedo_spk);
    if isempty(mrFet12_), continue; end
    vrDc2_site(iSite) = compute_dc2_(mrFet12_, viiSpk12_ord_, n1_, n2_, P); % Compute DC in CPU
    fprintf('.');    
end
dc = sqrt(abs(quantile(vrDc2_site, .5)));
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [vi, viSort] = rankorder_(vr, vcOrder)
% warning: 32 bit addressing
if nargin<2, vcOrder = 'ascend'; end
n=numel(vr);
[~,viSort] = sort(vr, vcOrder);
if isGpu_(vr)
    vi = zeros(n,1,'int32', 'gpuArray');
    vi(viSort) = 1:n;
else
    vi=zeros(n,1,'int32');
    vi(viSort) = 1:n;
end
end %func


%--------------------------------------------------------------------------
function plot_cdf_(vrSnr, fNorm)
if nargin<2, fNorm=0; end
vrX = sort(vrSnr,'ascend');
vrY = 1:numel(vrSnr);
if fNorm, vrY=vrY/vrY(end); end
stairs(vrY, vrX); 
end %func


%--------------------------------------------------------------------------
function [mr, miRange] = vr2mr3_(vr, vi, spkLim)
% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype

% prepare indices
if size(vi,2)==1, vi=vi'; end %row
viSpk = int32(spkLim(1):spkLim(end))';
miRange = bsxfun(@plus, viSpk, int32(vi));
miRange(miRange<1) = 1; 
miRange(miRange > numel(vr)) = numel(vr); %keep # sites consistent
% miRange = int32(miRange);

% build spike table
nSpks = numel(vi);
mr = reshape(vr(miRange(:)), [], nSpks);
end %func


%--------------------------------------------------------------------------
function [miSites, mrDist] = findNearSites_(mrSiteXY, maxSite, viSiteZero, viShank_site)
% find nearest sites
if nargin < 3, viSiteZero = []; end
if nargin<4, viShank_site = []; end
if numel(unique(viShank_site)) <= 1, viShank_site = []; end
if ~isempty(viSiteZero)
    mrSiteXY(viSiteZero,:) = inf; %bad sites will never be near
end
nNearSites = maxSite*2+1;
nSites = size(mrSiteXY,1);
nNearSites = min(nNearSites, nSites);
[miSites, mrDist] = deal(zeros(nNearSites, nSites));
viNearSites = 1:nNearSites;
for iSite=1:nSites
    vrSiteDist = pdist2_(mrSiteXY(iSite,:), mrSiteXY);
    if ~isempty(viShank_site)
        vrSiteDist(viShank_site(iSite) ~= viShank_site) = inf;
    end
    [vrSiteDist, viSrt] = sort(vrSiteDist, 'ascend');
    miSites(:,iSite) = viSrt(viNearSites);
    mrDist(:,iSite) = vrSiteDist(viNearSites);
end
end %func


%--------------------------------------------------------------------------
function vrDist_k = knn_sorted_(mrFet_srt, n_neigh, k_nearest)
nSpk = size(mrFet_srt, 1);
% if nargin<3, n_neigh=[]; end
% if nargin<4, k_nearest=[]; end

vrDist_k = zeros(nSpk, 1, 'single');
mrFet_srt = (mrFet_srt);
n1 = n_neigh*2+1;
mrFet_srt1 = mrFet_srt(1:n1,:);
iCirc = 1;
for iSpk = 1:nSpk
%     iSpk = viSpk_sub(iSpk1);
    if iSpk > n_neigh && iSpk <= nSpk-n_neigh
        mrFet_srt1(iCirc,:) = mrFet_srt(iSpk+n_neigh,:);
        iCirc=iCirc+1;
        if iCirc>n1, iCirc=1; end
    end
    vrDist1 = sort(sum(bsxfun(@minus, mrFet_srt1, mrFet_srt(iSpk,:)).^2, 2));
    vrDist_k(iSpk) = vrDist1(k_nearest);
end

vrDist_k = gather_(vrDist_k);
end %func


%--------------------------------------------------------------------------
function boxplot_(vrY, vrX, xbin, xlim1)
% range to plot: xlim1
vcMode = 'both';

viX = ceil(vrX/xbin);
ilim = ceil(xlim1/xbin);
viX(viX<ilim(1))=ilim(1);
viX(viX>ilim(end))=ilim(end);

nbins = diff(ilim)+1;
mrYp = zeros(nbins,3);
viXp = (ilim(1):ilim(end));
vrXp = viXp * xbin - xbin/2;
for ibin=1:nbins
    try
        ibin1 = viXp(ibin);
        mrYp(ibin, :) = quantile(vrY(viX==ibin1), [.25, .5, .75]);
    catch
        ;
    end
end

switch lower(vcMode)
    case 'stairs'
        vrXp = viXp - xbin/2;
        vrXp(end+1)=vrXp(end)+xbin;
        mrYp(end+1,:) = mrYp(end,:);
        stairs(vrXp, mrYp); grid on; 
    case 'line'
        plot(vrXp, mrYp); grid on; 
    case 'both'
        hold on; 
        vrXp1 = vrXp - xbin/2;
        vrXp1(end+1)=vrXp1(end)+xbin;      
        stairs(vrXp1, mrYp([1:end,end],[1, 3]), 'k-');
        
        plot(vrXp, mrYp(:,2), 'k.-', 'LineWidth', 1); 
        grid on; 
end
end %func


%--------------------------------------------------------------------------
function close_(hMsg)
try close(hMsg); catch; end
end


%--------------------------------------------------------------------------
function viTime1 = randomSelect_(viTime1, nShow)
if isempty(viTime1), return; end
if numel(viTime1) > nShow
    viTime1 = viTime1(randperm(numel(viTime1), nShow));
end
end %func


%--------------------------------------------------------------------------
function vr = linmap_(vr, lim1, lim2, fSat)
if nargin< 4
    fSat = 0;
end
if numel(lim1) == 1, lim1 = [-abs(lim1), abs(lim1)]; end
if numel(lim2) == 1, lim2 = [-abs(lim2), abs(lim2)]; end

if fSat
    vr(vr>lim1(2)) = lim1(2);
    vr(vr<lim1(1)) = lim1(1);
end
if lim1(1)==lim1(2)
    vr = vr / lim1(1);
else
    vr = interp1(lim1, lim2, vr, 'linear', 'extrap');
end
end %func


%--------------------------------------------------------------------------
function hPlot = plotTable_(lim, varargin)

vrX = floor((lim(1)*2:lim(2)*2+1)/2);
vrY = repmat([lim(1), lim(2), lim(2), lim(1)], [1, ceil(numel(vrX)/4)]);
vrY = vrY(1:numel(vrX));
hPlot = plot([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
end %func


%--------------------------------------------------------------------------
function hPlot = plotDiag_(lim, varargin)
[vrX, vrY] = plotDiag__(lim);
% vrY = floor((lim(1)*2:lim(2)*2+1)/2);
% vrX = [vrY(2:end), lim(end)];
% hPlot = plot([vrX(1:end-1), fliplr(vrY)], [vrY(1:end-1), fliplr(vrX)], varargin{:});
hPlot = plot(vrX, vrY, varargin{:});
end %func


%--------------------------------------------------------------------------
function [vrX, vrY] = plotDiag__(lim)
% lim: [start, end] or [start, end, offset]

vrY0 = floor((lim(1)*2:lim(2)*2+1)/2);
% vrY0 = lim(1):lim(end);
vrX0 = [vrY0(2:end), lim(2)];
vrX = [vrX0(1:end-1), fliplr(vrY0)];
vrY = [vrY0(1:end-1), fliplr(vrX0)];
if numel(lim)>=3
    vrX = vrX + lim(3);
    vrY = vrY + lim(3);
end
end %func


%--------------------------------------------------------------------------
function vlVisible = toggleVisible_(vhPlot, fVisible)
if isempty(vhPlot), return; end

if iscell(vhPlot)
    cvhPlot = vhPlot;
    if nargin<2
        cellfun(@(vhPlot)toggleVisible_(vhPlot), cvhPlot);
    else
        cellfun(@(vhPlot)toggleVisible_(vhPlot, fVisible), cvhPlot);
    end
    return;
end
try
    if nargin==1
        vlVisible = false(size(vhPlot));
        % toggle visibility
        for iH=1:numel(vhPlot)
            hPlot1 = vhPlot(iH);
            if strcmpi(get(hPlot1, 'Visible'), 'on')
                vlVisible(iH) = 0;
                set(hPlot1, 'Visible', 'off');
            else
                vlVisible(iH) = 1;
                set(hPlot1, 'Visible', 'on');
            end
        end
    else
        % set visible directly
        if fVisible
            vlVisible = true(size(vhPlot));
            set(vhPlot, 'Visible', 'on');
        else
            vlVisible = false(size(vhPlot));
            set(vhPlot, 'Visible', 'off');
        end
    end
catch
    return;
end
end %func


%--------------------------------------------------------------------------
function S = struct_delete_(S, varargin)
% delete and set to empty

for i=1:numel(varargin)
    try 
        delete(S.(varargin{i}));
        S.(varargin{i}) = [];
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function hMsg = msgbox_open_(vcMessage)
global fDebug_ui

% if get_set_([], 'fDebug_ui', 0), hMsg = []; return; end
if fDebug_ui==1, hMsg = []; return; end
hFig = gcf;
hMsg = msgbox(vcMessage);  
figure(hFig); drawnow;
end


%--------------------------------------------------------------------------
% function S = struct_set_(S, index, val, varargin)
% for i=1:numel(varargin)
%     vcVar = varargin{i};
%     if ~isfield(S, vcVar), continue; end %ignore if not
%     var1 = S.(vcVar);
%     if isempty(val)
%         var1(index) = [];
%     else
%         if iscell(var1)
%             var1{index} = val;
%         else
%             var1(index) = val;
%         end
%     end
%     S.(vcVar) = var1;
% end
% end %func


%--------------------------------------------------------------------------
function [i1, i2] = swap_(i1, i2)
i11 = i1;   i1 = i2;    i2 = i11;
end %func


%--------------------------------------------------------------------------
function Sclu = merge_clu_pair_(Sclu, iClu1, iClu2)
% if iClu1>iClu2, [iClu1, iClu2] = swap(iClu1, iClu2); end

% update vnSpk_clu, viClu, viSite_clu. move iClu2 to iClu1
n1 = Sclu.vnSpk_clu(iClu1);
n2 = Sclu.vnSpk_clu(iClu2);
Sclu.vnSpk_clu(iClu1) = n1 + n2;
Sclu.vnSpk_clu(iClu2) = 0;
Sclu.viClu(Sclu.viClu == iClu2) = iClu1;
Sclu.cviSpk_clu{iClu1} = find(Sclu.viClu == iClu1);
Sclu.cviSpk_clu{iClu2} = [];
try
    Sclu.csNote_clu{iClu1} = '';
    Sclu.csNote_clu{iClu2} = '';
catch
end
end %func


%--------------------------------------------------------------------------
function set_axis_(hFig, xlim1, ylim1, xlim0, ylim0)
% set the window within the box limit
if nargin <= 3
    % square case
    xlim0 = ylim1;
    ylim1 = xlim1;
    ylim0 = xlim0;
end

hFig_prev = gcf;
figure(hFig); 
dx = diff(xlim1);
dy = diff(ylim1);

if xlim1(1)<xlim0(1), xlim1 = xlim0(1) + [0, dx]; end
if xlim1(2)>xlim0(2), xlim1 = xlim0(2) + [-dx, 0]; end
if ylim1(1)<ylim0(1), ylim1 = ylim0(1) + [0, dy]; end
if ylim1(2)>ylim0(2), ylim1 = ylim0(2) + [-dy, 0]; end

xlim1(1) = max(xlim1(1), xlim0(1));
ylim1(1) = max(ylim1(1), ylim0(1));
xlim1(2) = min(xlim1(2), xlim0(2));
ylim1(2) = min(ylim1(2), ylim0(2));

axis([xlim1, ylim1]);

figure(hFig_prev);
end %func


%--------------------------------------------------------------------------
function mrDist_clu = S_clu_wavcor_1_(S_clu, P)
% oldd version
% viShift0 = 0; %-2:2; %compare time shifted version
viShift0 = -6:6; %time shift range
fWav_cumsum = 0;

% tmrWav_clu = ifeq_(fWav_cumsum, S_clu.tmrWav_clu, S_clu.tmrWav_spk_clu);
tmrWav_clu = meanSubt_(S_clu.tmrWav_raw_clu);
% maxSite = P.maxSite_merge;
maxSite = ceil(P.maxSite/2);
mrDist_clu = nan(S_clu.nClu);
miSites = P.miSites;
% miSites = P.miSites(1:(maxSite+1), :); % center
for iClu2 = 1:S_clu.nClu         
    iSite2 = S_clu.viSite_clu(iClu2);   
    viSite2 = miSites(:,iSite2);
    vrWav2 = get_wav_(tmrWav_clu, viSite2, iClu2);
    mrWav2 = zscore_(shift_vr_(vrWav2(:), viShift0)); %try +/-4 samples
    viClu1 = find(abs(S_clu.viSite_clu - iSite2) <= maxSite);
    viClu1(viClu1 == iClu2) = [];
    viClu1 = viClu1(:)';
    vlWav2 = vrWav2~=0;
    for iClu1 = viClu1   
        try            
            vrWav1 = zscore_(get_wav_(tmrWav_clu, viSite2, iClu1));
            n12 = sum(vlWav2 & vrWav1~=0);
            corr12 = max(vrWav1' * mrWav2) / n12; %numel(vrWav1);
            mrDist_clu(iClu1, iClu2)=corr12;
        catch
            disp('Sclu_wavcor: error');
        end
    end
end
end %func


%--------------------------------------------------------------------------
function mr = get_wav_(tmrWav_clu, viSite, iClu, viShift)
if nargin <4, viShift = 0; end

mrWav1 = (tmrWav_clu(:,viSite, iClu));
% mrWav1 = zscore_(mrWav1);

% induce waveform shift
cmr = cell(1, numel(viShift));
for iShift = 1:numel(viShift)
	mr1 = shift_mr_(mrWav1, viShift(iShift));    
    mr1 = bsxfun(@minus, mr1, mean(mr1,1));
    cmr{iShift} = mr1(:);
end
% mr = zscore_(cell2mat_(cmr));
mr = cell2mat_(cmr);
end %func


%--------------------------------------------------------------------------
function mr = shift_mr_(mr, n)
if n<0
    n=-n;
    mr(1:end-n, :) = mr(n+1:end, :);
elseif n>0
    mr(n+1:end, :) = mr(1:end-n, :);    
else %n==0
    return;     
end
end %func


%--------------------------------------------------------------------------
% function [S_clu, nMerged] = S_clu_merge_mtx_(S_clu, mrDist, dist_thresh)
% % merge less than
% % merge clusters based on the distance creterion
% nMerged = 0;
% fMutualMin = 0; % setting this reduces false positive rate
% if isempty(dist_thresh), return; end
% fprintf('S_clu_merge_mtx: merging based on waveform similarity.\n\t');
% t1=tic;
% n = size(mrDist,2);
% viClu_merge = 1:n;
% vlKeep_clu = true(n, 1);
% if fMutualMin
%     [inds1,inds2] = get_pairs_to_compare_(mrDist, dist_thresh);
%     viClu_merge(inds1) = inds2;
%     vlKeep_clu(inds2) = 0;
% else
%     %mrDist = (mrDist + mrDist') / 2;
%     mrDist = max(mrDist, mrDist');
%     mrDist(tril(true(n))) = nan; %ignore bottom half
%     % mrDist(mrDist==0) = nan;        
%     for iClu = 1:size(mrDist,2)
%         [dist1, iClu_merge] = min(mrDist(:,iClu));
%         if dist1 <= dist_thresh
%             viClu_merge(iClu) = viClu_merge(iClu_merge);
%             vlKeep_clu(iClu) = 0;
%         end
%         fprintf('.');
%     end
% end
% 
% % Merge clusters
% [~, vi_c2b, vi_b2c] = unique(viClu_merge);
% % reassign cluster number
% nClu_old = S_clu.nClu;
% vlPos = S_clu.viClu>0;
% S_clu.viClu(vlPos) = vi_b2c(viClu_merge(S_clu.viClu(vlPos)));
% try
% S_clu.icl(~vlKeep_clu) = [];
% catch
%     ;
% end
% S_clu.nClu = max(vi_b2c);
% S_clu = S_clu_refresh_(S_clu);
% fprintf('\n\tnClu: %d->%d, took %0.1fs\n', nClu_old, S_clu.nClu, toc(t1));
% nMerged = nClu_old - S_clu.nClu;
% end %function


%--------------------------------------------------------------------------
function [S_clu, nClu_merged] = S_clu_wavcor_merge_(S_clu, P)
global tnWav_spk tnWav_raw

[viSite_spk] = get0_('viSite_spk');
mrWavCor = S_clu.mrWavCor;
nClu = size(mrWavCor, 2);

% Identify clusters to remove, update and same (no change), disjoint sets
% mrWavCor(tril(true(nClu)) | mrWavCor==0) = nan; %ignore bottom half
mrWavCor(tril(true(nClu))) = 0; %ignore bottom half
[vrCor_max, viMap_clu] = max(mrWavCor);
vlKeep_clu = vrCor_max < P.maxWavCor; % | isnan(vrCor_max);
if all(vlKeep_clu), nClu_merged=0; return ;end
min_cor = min(vrCor_max(~vlKeep_clu));
viClu_same = find(vlKeep_clu);
viMap_clu(vlKeep_clu) = viClu_same;
viClu_same = setdiff(viClu_same, viMap_clu(~vlKeep_clu));
viClu_remove = setdiff(1:nClu, viMap_clu);
viClu_update = setdiff(setdiff(1:nClu, viClu_same), viClu_remove);
% viClu_update = setdiff(1:nClu, viClu_same);

% update cluster number
try S_clu.icl(viClu_remove) = []; catch, end
S_clu = S_clu_map_index_(S_clu, viMap_clu); %index mapped
P.fVerbose = 0;
S_clu = S_clu_refrac_(S_clu, P); % remove refrac spikes

% update cluster waveforms and distance
S_clu = clu2wav_(S_clu, viSite_spk, tnWav_spk, tnWav_raw, viClu_update); %update cluster waveforms
S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update);
% S_clu = S_clu_refresh_(S_clu); % remove empty and remap
S_clu = S_clu_remove_empty_(S_clu);

nClu_merged = nClu - S_clu.nClu;
fprintf('\n\tnClu: %d->%d (%d merged, min-cor: %0.4f)\n', nClu, S_clu.nClu, nClu_merged, min_cor);
end %func


%--------------------------------------------------------------------------
function varargout = multiindex_(viKeep, varargin)
% index first dimension of variables by a given index

if nargout ~= numel(varargin), error('Number of argin=argout'); end
    
if islogical(viKeep), viKeep = find(viKeep); end
for i=1:numel(varargin)
    var1 = varargin{i};
    if isvector(var1)
        varargout{i} = var1(viKeep);
    elseif ismatrix(var1)
        varargout{i} = var1(viKeep,:);
    else %multidimensional variable
        varargout{i} = var1(viKeep,:,:);
    end    
end
end %func


%--------------------------------------------------------------------------
function S = struct_reorder_(S, viKeep, varargin)
for i=1:numel(varargin)
    try
        vcVar = varargin{i};
        if ~isfield(S, vcVar), continue; end %ignore if not
        vr1 = S.(vcVar);
        if isvector(vr1)
            vr1 = vr1(viKeep);
        elseif ismatrix(vr1)
            vr1 = vr1(viKeep, :);
        else
            vr1 = vr1(viKeep, :, :);
        end
        S.(vcVar) = vr1;
    catch
        ;
    end
end
end %func


%--------------------------------------------------------------------------
function [ vlSpkIn ] = auto_split_wav_(mrSpkWav)
% TODO: ask users number of clusters and split multi-way
%Make automatic split of clusters using PCA + kmeans clustering
%  input Sclu, trSpkWav and cluster_id of the cluster you want to cut

% show clusters and ask number of split
nSplits = 2;
nSpks = size(mrSpkWav,2);
[~,mrFet,~] = pca(double(mrSpkWav'), 'Centered', 1, 'NumComponents', 3);

% nSplit = preview_split_(mrSpkWav1);
% if isnan(nSplit), return; end

try    
    % kmean clustering into 2
    idx = kmeans(mrFet, nSplits);     
    mrFet = mrFet';
    d12 = mad_dist_(mrFet(:,idx==1), mrFet(:,idx==2)); 
    fprintf('mad_dist: %f\n', d12);
%     idx = kmeans([pca_1,pca_2], NUM_SPLIT);
    vlSpkIn = logical(idx-1);
catch
%         msgbox('Too few spikes to auto-split');
    vlSpkIn = false(nSpks,1);
    vlSpkIn(1:end/2) = true;
    return;
end



resize_figure_([], [.5 0 .5 1]);

subplot(2,3,1:3)
hold on
plot(mrFet(1,vlSpkIn),mrFet(2,vlSpkIn),'b.')
plot(mrFet(1,~vlSpkIn),mrFet(2,~vlSpkIn),'r.')
xlabel('PCA 1st Component')
ylabel('PCA 2nd Component')
hold off

min_y=min(reshape(mrSpkWav,1,[]));
max_y=max(reshape(mrSpkWav,1,[]));
subplot(2,3,4)
hold on
title('Spikes of Cluster 1')
viSpk1 = subsample_vr_(1:size(mrSpkWav,1), 1000); % show 1000 sites only
plot(mrSpkWav(viSpk1,vlSpkIn))
ylim([min_y max_y])
hold off

subplot(2,3,5)
hold on
title('Spikes of Cluster 2')
plot(mrSpkWav(viSpk1,~vlSpkIn))
ylim([min_y max_y])
hold off

subplot(2,3,6)
hold on
title('Overlay of mean spike')
plot(mean(mrSpkWav(viSpk1,vlSpkIn),2),'b')
plot(mean(mrSpkWav(viSpk1,~vlSpkIn),2),'r')
ylim([min_y max_y])
hold off
end %func


%--------------------------------------------------------------------------
function [vi, nClu, viA] = mapIndex_(vi, viA, viB)
% change the index of vi according to the map (viA)

if nargin<2, viA = setdiff(unique(vi), 0); end %excl zero
if nargin<3, viB = 1:numel(viA); end
nClu = viB(end);
viAB(viA) = viB; %create a translation table A->B
vl = vi>0;
vi(vl) = viAB(vi(vl)); %do not map zeros
end %func


%--------------------------------------------------------------------------
function d12 = mad_dist_(mrFet1, mrFet2)
% distance between two clusters
if ~ismatrix(mrFet1)
    mrFet1 = reshape(mrFet1, [], size(mrFet1,3));
end
if ~ismatrix(mrFet2)
    mrFet2 = reshape(mrFet2, [], size(mrFet2,3));
end
vrFet1_med = median(mrFet1, 2);
vrFet2_med = median(mrFet2, 2);
vrFet12_med = vrFet1_med - vrFet2_med;
norm12 = sqrt(sum(vrFet12_med.^2));
vrFet12_med1 = vrFet12_med / norm12;
mad1 = median(abs(vrFet12_med1' * bsxfun(@minus, mrFet1, vrFet1_med)));
mad2 = median(abs(vrFet12_med1' * bsxfun(@minus, mrFet2, vrFet2_med)));
d12 = norm12 / sqrt(mad1.^2 + mad2.^2);
end %func


%--------------------------------------------------------------------------
function mr = filterq_(vrA, vrB, mr, dimm)
% quick filter using single instead of filtfilt
% faster than filtfilt and takes care of the time shift 

if nargin < 4, dimm = 1; end
if numel(vrA)==1, return; end
if isempty(vrB), vrB=sum(vrA); end
%JJJ 2015 09 16
% mr = filter(vrA, vrB, mr, [], dimm);
mr = circshift(filter(vrA, vrB, mr, [], dimm), -ceil(numel(vrA)/2), dimm);
end %func


%--------------------------------------------------------------------------
function S = imec3_imroTbl_(cSmeta)
% Smeta has imroTbl

vcDir_probe = 'C:\Dropbox (HHMI)\IMEC\SpikeGLX_Probe_Cal_Data\';  %this may crash. probe calibaration folder

if isstruct(cSmeta), cSmeta = {cSmeta}; end %turn it into cell of struct
% parse imroTbl
cs_imroTbl = cellfun(@(S)S.imroTbl, cSmeta, 'UniformOutput', 0);
cvn_imroTbl = cellfun(@(vc)textscan(vc, '%d', 'Delimiter', '( ),'), cs_imroTbl, 'UniformOutput', 0);
cvn_imroTbl = cellfun(@(c)c{1}, cvn_imroTbl, 'UniformOutput', 0);   
S.viBank = cellfun(@(vn)vn(7), cvn_imroTbl);
S.viRef = cellfun(@(vn)vn(8), cvn_imroTbl);
S.vrGain_ap = single(cellfun(@(vn)vn(9), cvn_imroTbl));
S.vrGain_lf = single(cellfun(@(vn)vn(10), cvn_imroTbl));
S.nSites_bank = cvn_imroTbl{1}(4);

Smeta1 = cSmeta{1};

% correct gain
nFiles = numel(S.viBank);
nSites = numel(Smeta1.viSites);
[mrScale_ap, mrScale_lf] = deal(ones(nSites, nFiles));
S.vcProbeSN = sprintf('1%d%d', Smeta1.imProbeSN, Smeta1.imProbeOpt);
% read gain correction
vrGainCorr = ones(1, S.nSites_bank*4);
if Smeta1.imProbeOpt ~= 2
    try    
        vcFile_csv = sprintf('%s1%d%d\\Gain correction.csv', vcDir_probe, Smeta1.imProbeSN, Smeta1.imProbeOpt);
        try
            vrGainCorr = csvread(vcFile_csv, 1, 1);
        catch
            vrGainCorr = csvread(vcFile_csv, 0, 0);
        end
    catch
        ;
    end
end

% build scale
for iFile = 1:nFiles
    vrGainCorr1 = vrGainCorr(Smeta1.viSites + double(S.nSites_bank*S.viBank(iFile)));
    mrScale_ap(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_ap(iFile) .* vrGainCorr1;
    mrScale_lf(:,iFile) = 1.2 * 1e6 / 2^10 / S.vrGain_lf(iFile) .* vrGainCorr1;
end
S.mrScale_ap = mrScale_ap;
S.mrScale_lf = mrScale_lf;
end %func


%--------------------------------------------------------------------------
function plot_raster_(P, iClu, S_clu)
%plot_raster_()
%   plot if window open using curretnly selected clusters
%plot_raster_(P, iClu, S_clu)
%   Open window and plot specific clusters and S_clu

persistent hFig

% import  trial time
% P = loadParam(vcFile_prm);
if nargin==0
    if ~isvalid_(hFig), return; end
    S0 = get0_();
    P = S0.P; iClu = [S0.iCluCopy, S0.iCluPaste]; S_clu = S0.S_clu;
else
    if nargin < 2, iClu = []; end
    if nargin < 3, S_clu = []; end
end

if isempty(S_clu), S_clu = get0_('S_clu'); end

if isfield(P, 'vcFile_psth'), P.vcFile_trial = P.vcFile_psth; end
%
try
    crTime_trial = loadTrial_(P.vcFile_trial);
catch
    return;
end
%viTime_spk = get0_('viTime_spk');
%vrTime_trial = loadTrial_(P.vcFile_trial);
if ~iscell(crTime_trial)
    crTime_trial = {crTime_trial};
end
nstims = numel(crTime_trial);
if isempty(crTime_trial), msgbox('Trial file does not exist', 'modal'); return; end
if ~isvalid_(hFig)
    hFig = create_figure_('FigTrial', [.5  0 .5 1], P.vcFile_trial, 0, 0);
end
%figure(hFig); clf(hFig); 
hTabGroup = uitabgroup(hFig);
% offset = 0;
if isempty(iClu)
    viClu_plot = 1:S_clu.nClu;
else
    viClu_plot = iClu; %copy and paste
end
for iClu = viClu_plot
    htab1 = uitab(hTabGroup, 'Title', sprintf('Clu %d', iClu));    
    axoffset = 0.05;
    axlen = 0.9/nstims;
    [vhAx1, vhAx2] = deal(nan(nstims, 1));
    for iStim = 1:nstims
        vrTime_trial = crTime_trial{iStim}; %(:,1);
        nTrials = numel(vrTime_trial);
        
        % plot raster
        vhAx1(iStim) = axes('Parent', htab1, 'Position',[.08 axoffset .9 axlen*.68]);
        viTime_clu1 = S_clu_time_(S_clu, iClu);
        plot_raster_clu_(viTime_clu1, vrTime_trial, P, vhAx1(iStim));
        
        % plot psth
        vhAx2(iStim) = axes('Parent', htab1, 'Position',[.08 axoffset + axlen*.68 .9 axlen*.2]);
        plot_psth_clu_(viTime_clu1, vrTime_trial, P, vhAx2(iStim));
        axoffset = axoffset + axlen;
        
        title(vhAx2(iStim), sprintf('Cluster %d; %d trials', iClu, nTrials));
    end
%     offset = offset + nTrials;
    if numel(vhAx1)>2
        set(vhAx1(2:end),'xticklabel',{});
        for ax = vhAx1(2:end)
            xlabel(ax, '')
        end
        
    end
end
end %func


%--------------------------------------------------------------------------
function plot_psth_clu_(viTime_clu, vrTime_trial, P, hAx)
if nargin<4, hAx=gca; end
tbin = P.tbin_psth;
nbin = round(tbin * P.sRateHz);
nlim = round(P.tlim_psth/tbin);
viTime_Trial = round(vrTime_trial / tbin);

vlTime1=zeros(0);
vlTime1(ceil(double(viTime_clu)/nbin))=1;
mr1 = vr2mr2_(double(vlTime1), viTime_Trial, nlim);
vnRate = mean(mr1,2) / tbin;
vrTimePlot = (nlim(1):nlim(end))*tbin + tbin/2;

bar(hAx, vrTimePlot, vnRate, 1, 'EdgeColor', 'none');

vrXTick = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
set(hAx, 'XTick', vrXTick, 'XTickLabel', []);
grid(hAx, 'on');
hold(hAx, 'on'); 
plot(hAx, [0 0], get(hAx,'YLim'), 'r-');
ylabel(hAx, 'Rate (Hz)');
xlim(hAx, P.tlim_psth);
end


%--------------------------------------------------------------------------
function plot_raster_clu_(viTime_clu, vrTime_trial, P, hAx)
if nargin<4, hAx=gca; end

trialLength = diff(P.tlim_psth); % seconds
nTrials = numel(vrTime_trial);
spikeTimes = cell(nTrials, 1);
t0 = -P.tlim_psth(1);
for iTrial = 1:nTrials
    rTime_trial1 = vrTime_trial(iTrial);
    vrTime_lim1 = rTime_trial1 + P.tlim_psth;
    vrTime_clu1 = double(viTime_clu) / P.sRateHz;
    vrTime_clu1 = vrTime_clu1(vrTime_clu1>=vrTime_lim1(1) & vrTime_clu1<vrTime_lim1(2));
    
    vrTime_clu1 = (vrTime_clu1 - rTime_trial1 + t0) / trialLength;
    spikeTimes{iTrial} = vrTime_clu1';
end

% Plot
hAx_pre = axes_(hAx);
plotSpikeRaster(spikeTimes,'PlotType','vertline','RelSpikeStartTime',0.01,'XLimForCell',[0 1], ...
    'LineFormat', struct('LineWidth', 1.5));
axes_(hAx_pre);
ylabel(hAx, 'Trial #')
% title('Vertical Lines With Spike Offset of 10ms (Not Typical; for Demo Purposes)');
vrXTickLabel = P.tlim_psth(1):(P.xtick_psth):P.tlim_psth(2);
vrXTick = linspace(0,1,numel(vrXTickLabel));
set(hAx, {'XTick', 'XTickLabel'}, {vrXTick, vrXTickLabel});
grid(hAx, 'on');
hold(hAx, 'on'); 
plot(hAx, [t0,t0]/trialLength, get(hAx,'YLim'), 'r-');
xlabel(hAx, 'Time (s)');
end %func


%--------------------------------------------------------------------------
function mr = vr2mr2_(vr, viRow, spkLim, viCol)
if nargin<4, viCol = []; end
% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype

% prepare indices
if size(viRow,2)==1, viRow=viRow'; end %row
viSpk = int32(spkLim(1):spkLim(end))';
miRange = bsxfun(@plus, viSpk, int32(viRow));
miRange = min(max(miRange, 1), numel(vr));
if isempty(viCol)
    mr = vr(miRange); %2x faster
else
    mr = vr(miRange, viCol); %matrix passed to save time
end
end %func


%--------------------------------------------------------------------------
function hAx_prev = axes_(hAx)
hAx_prev = gca;
axes(hAx);
end


%--------------------------------------------------------------------------
function [viTime_spk, vrAmp_spk, viSite_spk] = spike_refrac_(viTime_spk, vrAmp_spk, viSite_spk, nRefrac)
% Remove smaller spikes if a bigger one detected within nRefrac
% spike_refrac_(viSpk, vrSpk, [], nRefrac)
% spike_refrac_(viSpk, vrSpk, viSite, nRefrac)
nRepeat_max = 10;
for iRepeat = 1:nRepeat_max
    if isempty(viTime_spk), return ;end
    vnDiff_spk = diff(diff(viTime_spk) <= nRefrac);
    viStart = find(vnDiff_spk > 0) + 1;
    viEnd = find(vnDiff_spk < 0) + 1;
    if isempty(viStart) || isempty(viEnd), return; end
    viEnd(viEnd < viStart(1)) = [];
    if isempty(viEnd), return; end
    viStart(viStart > viEnd(end)) = [];
    nGroup = numel(viStart);
    if nGroup==0, return; end
    assert(nGroup == numel(viEnd), 'spike_refrac_:nStart==nEnd');    
    vlRemove = false(size(viTime_spk));
    
    for iGroup = 1:nGroup
        vi_ = viStart(iGroup):viEnd(iGroup);
        [~, imax_] = max(vrAmp_spk(vi_));
        vlRemove(vi_(imax_)) = true;
    end
    viTime_spk(vlRemove) = [];
    if ~isempty(vrAmp_spk), vrAmp_spk(vlRemove) = []; end
    if ~isempty(viSite_spk), viSite_spk(vlRemove) = []; end
end
end %func


%--------------------------------------------------------------------------
function [viSpk, vrSpk, viSite] = spike_refrac__(viSpk, vrSpk, viSite, nRefrac)
% Remove smaller spikes if a bigger one detected within nRefrac
% spike_refrac_(viSpk, vrSpk, [], nRefrac)
% spike_refrac_(viSpk, vrSpk, viSite, nRefrac)

nSkip_refrac = 8;
% remove refractory period
vlKeep = true(size(viSpk));
% if isGpu_(viSpk), vlKeep = gpuArray(vlKeep); end
while (1)
    viKeep1 = find(vlKeep);
    viRefrac1 = find(diff(viSpk(viKeep1)) <= nRefrac);
    if isempty(viRefrac1), break; end
    
    vi1 = viRefrac1(1:nSkip_refrac:end);     
    viRemoveA = viKeep1(vi1);
    viRemoveB = viKeep1(vi1+1);   
    if ~isempty(vrSpk)
        vl1 = abs(vrSpk(viRemoveA)) < abs(vrSpk(viRemoveB));
        vlKeep(viRemoveA(vl1)) = 0;
        vlKeep(viRemoveB(~vl1)) = 0;
    else
        vlKeep(viRemoveB) = 0;
    end
end

viSpk(~vlKeep) = [];
if ~isempty(vrSpk), vrSpk(~vlKeep) = []; end
if ~isempty(viSite), viSite(~vlKeep) = []; end
end %func


%--------------------------------------------------------------------------
function export_jrc1_(vcFile_prm)
% Export to version 1 format (_clu.mat) and (_evt.mat)
% error('export_jrc1_: not implemented yet');
% Load info from previous version: time, site, spike
P = loadParam_(vcFile_prm);
S0 = load(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));

fprintf('Exporting to JRCLUST ver.1 format\n');
% Build Sevt and save
Sevt = S0;
Sevt = rmfield_(Sevt, 'S_clu', 'cmrFet_site', 'cvrTime_site', 'cvrVpp_site', ...
    'dimm_raw', 'dimm_spk', 'miSites_fet', 'mrFet', 'runtime_detect', 'runtime_sort', ...
    'viSite_spk', 'viT_offset_file', 'viTime_spk', 'vrAmp_spk');
Sevt.cvrSpk_site = vr2cell_(S0.vrAmp_spk, S0.cviSpk_site); %S0.cvrVpp_site;
% Sevt.miSites_fet = S0.miSites_fet;
[Sevt.mrPv, Sevt.mrWav_spk, Sevt.trFet] = deal([]);
Sevt.viSite = S0.viSite_spk;
Sevt.viSpk = S0.viTime_spk;
Sevt.vrSpk = S0.vrAmp_spk;
Sevt.vrThresh_uV = bit2uV_(S0.vrThresh_site, S0.P);
Sevt.dimm_fet = [1, S0.dimm_fet(:)'];
write_struct_(strrep(P.vcFile_prm, '.prm', '_evt.mat'), Sevt);
fprintf('\tEvent struct (Sevt) exported.\n\t');
assignWorkspace_(Sevt);

% S_clu
if isfield(S0, 'S_clu')
    Sclu = S0.S_clu;
    Sclu = rmfield_(Sclu, 'tmrWav_raw_clu', 'tmrWav_spk_clu', 'trWav_raw_clu', 'trWav_spk_clu', 'viSite_min_clu', 'vrVmin_clu');
    Sclu.trWav_dim = [];
    Sclu.viSite = S0.viSite_spk;
    Sclu.viTime = S0.viTime_spk;
    Sclu.vrSnr_Vmin_clu = S0.S_clu.vrSnr_clu;
    write_struct_(strrep(P.vcFile_prm, '.prm', '_clu.mat'), Sclu);
    fprintf('\tCluster struct (Sclu) exported.\n\t');
    assignWorkspace_(Sclu);
end
end %func


%--------------------------------------------------------------------------
function cvr = vr2cell_(vr, cvi)
cvr = cellfun(@(vi)vr(vi), cvi, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function [S_clu, S0] = S_clu_new_(arg1, S0)
% S_clu = S_clu_new_(S_clu, S0)
% S_clu = S_clu_new_(viClu, S0)

global tnWav_spk tnWav_raw
if nargin<2, S0 = get(0, 'UserData'); end
% S_clu = S0.S_clu;
S_clu = get_(S0, 'S_clu'); %previous S_clu
if isempty(S_clu), S_clu = struct([]); end
if ~isstruct(arg1)
    S_clu.viClu = arg1; %skip FigRD step for imported cluster
else
    S_clu = struct_append_(S_clu, arg1);
end
S_clu.viClu = int32(S_clu.viClu);
S_clu = S_clu_refresh_(S_clu);
S_clu = S_clu_update_wav_(S_clu, S0, S0.P);
S_clu = S_clu_position_(S_clu);
if ~isfield(S_clu, 'csNote_clu')
    S_clu.csNote_clu = cell(S_clu.nClu, 1);  %reset note
end
S_clu = S_clu_quality_(S_clu, S0.P);
[S_clu, S0] = S_clu_commit_(S_clu, 'S_clu_new_');
end %func


%--------------------------------------------------------------------------
function [mr, vi_shuffle] = shuffle_static_(mr, dimm)
% dimm = 1 or 2 (dimension to shuffle
if nargin<2, dimm=1; end
fStatic = 0;
if fStatic
    s = RandStream('mt19937ar','Seed',0); %always same shuffle order
    switch dimm
        case 1
            vi_shuffle = randperm(s, size(mr,1));
            mr = mr(vi_shuffle, :);
        case 2
            vi_shuffle = randperm(s, size(mr,2));
            mr = mr(:, vi_shuffle);
    end
else
    switch dimm
        case 1
            mr = mr(randperm(size(mr,1)), :);
        case 2
            mr = mr(:, randperm(size(mr,2)));
    end
end
end %func


%--------------------------------------------------------------------------
function [allScores, allFPs, allMisses, allMerges] = compareClustering2_(cluGT, resGT, cluTest, resTest)
% function compareClustering(cluGT, resGT, cluTest, resTest[, datFilename])
% - clu and res variables are length nSpikes, for ground truth (GT) and for
% the clustering to be evaluated (Test). 
% kilosort, Marius Pachitariu, 2016-dec-21

t1 = tic; fprintf('kilosort-style ground truth validation\n\t');
[cluGT, cluTest, resTest] = multifun_(@double, cluGT, cluTest, resTest);
resGT = int64(resGT);
GTcluIDs = unique(cluGT);
testCluIDs = unique(cluTest);
jitter = 12;

nSp = zeros(max(testCluIDs), 1);
for j = 1:max(testCluIDs)
    nSp(j) = max(1, sum(cluTest==j));
end
nSp0 = nSp;

for cGT = 1:length(GTcluIDs)
    rGT = int32(resGT(cluGT==GTcluIDs(cGT)));
    S = spalloc(numel(rGT), max(testCluIDs), numel(rGT) * 10);
    % find the initial best match
    mergeIDs = [];
    scores = [];
    falsePos = [];
    missRate = [];
    
    igt = 1;
    
    nSp = nSp0;
    nrGT = numel(rGT);
    flag = false;
    for j = 1:numel(cluTest)
        while (resTest(j) > rGT(igt) + jitter)
            % the curent spikes is now too large compared to GT, advance the GT
            igt = igt + 1;
            if igt>nrGT
               flag = true;
               break;
            end
        end
        if flag, break; end        
        if resTest(j)>rGT(igt)-jitter
            % we found a match, add a tick to the right cluster
              S(igt, cluTest(j)) = 1;
        end
    end
    numMatch = sum(S,1)';
    misses = (nrGT-numMatch)/nrGT; % missed these spikes, as a proportion of the total true spikes
    fps = (nSp-numMatch)./nSp; % number of comparison spikes not near a GT spike, as a proportion of the number of guesses

    sc = 1-(fps+misses);
    best = find(sc==max(sc),1);
    mergeIDs(end+1) = best;
    scores(end+1) = sc(best);
    falsePos(end+1) = fps(best);
    missRate(end+1) = misses(best);
    
%     fprintf(1, '  found initial best %d: score %.2f (%d spikes, %.2f FP, %.2f miss)\n', ...
%         mergeIDs(1), scores(1), sum(cluTest==mergeIDs(1)), fps(best), misses(best));
    
    S0 = S(:, best);
    nSp = nSp + nSp0(best);
    while scores(end)>0 && (length(scores)==1 || ( scores(end)>(scores(end-1) + 1*0.01) && scores(end)<=0.99 ))
        % find the best match
        S = bsxfun(@max, S, S0);
        
        numMatch = sum(S,1)';
        misses = (nrGT-numMatch)/nrGT; % missed these spikes, as a proportion of the total true spikes
        fps = (nSp-numMatch)./nSp; % number of comparison spikes not near a GT spike, as a proportion of the number of guesses
        
        sc = 1-(fps+misses);
        best = find(sc==max(sc),1);
        mergeIDs(end+1) = best;
        scores(end+1) = sc(best);
        falsePos(end+1) = fps(best);
        missRate(end+1) = misses(best);
        
%         fprintf(1, '    best merge with %d: score %.2f (%d/%d new/total spikes, %.2f FP, %.2f miss)\n', ...
%             mergeIDs(end), scores(end), nSp0(best), nSp(best), fps(best), misses(best));
        
        S0 = S(:, best);
        nSp = nSp + nSp0(best);
                
    end
    
    if length(scores)==1 || scores(end)>(scores(end-1)+0.01)
        % the last merge did help, so include it
        allMerges{cGT} = mergeIDs(1:end);
        allScores{cGT} = scores(1:end);
        allFPs{cGT} = falsePos(1:end);
        allMisses{cGT} = missRate(1:end);
    else
        % the last merge actually didn't help (or didn't help enough), so
        % exclude it
        allMerges{cGT} = mergeIDs(1:end-1);
        allScores{cGT} = scores(1:end-1);
        allFPs{cGT} = falsePos(1:end-1);
        allMisses{cGT} = missRate(1:end-1);
    end
    fprintf('.');
end

initScore = zeros(1, length(GTcluIDs));
finalScore = zeros(1, length(GTcluIDs));
numMerges = zeros(1, length(GTcluIDs));
fprintf(1, 'Validation score (Kilosort-style)\n'); 
for cGT = 1:length(GTcluIDs)
     initScore(cGT) = allScores{cGT}(1);
     finalScore(cGT) = allScores{cGT}(end);
     numMerges(cGT) = length(allScores{cGT})-1;
end
finalScores = cellfun(@(x) x(end), allScores);
nMerges = cellfun(@(x) numel(x)-1, allMerges);

fprintf(1, '\tmedian initial score: %.2f; median best score: %.2f\n', median(initScore), median(finalScore));
fprintf(1, '\ttotal merges required: %d\n', sum(numMerges));
fprintf('\t%d / %d good cells, score > 0.8 (pre-merge) \n', sum(cellfun(@(x) x(1), allScores)>.8), numel(allScores))
fprintf('\t%d / %d good cells, score > 0.8 (post-merge) \n', sum(cellfun(@(x) x(end), allScores)>.8), numel(allScores))
fprintf('\tMean merges per good cell %2.2f \n', mean(nMerges(finalScores>.8)))
fprintf('\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function S_ksort = kilosort_(vcFile_prm)
% Run Kilosort
fSavePhy = 1;
fMerge_post = 1;
if ischar(vcFile_prm)
    fprintf('Running kilosort on %s\n', vcFile_prm); 
    P = loadParam_(vcFile_prm);
    if isempty(P), return; end
else
    P = vcFile_prm;
    vcFile_prm = P.vcFile_prm;
    fprintf('Running kilosort on %s\n', vcFile_prm); 
end
runtime_ksort = tic;

% Run kilosort
[fpath, ~, ~] = fileparts(vcFile_prm);
ops = kilosort('config', P); %get config
S_chanMap = kilosort('chanMap', P); %make channel map

[rez, DATA, uproj] = kilosort('preprocessData', ops); % preprocess data and extract spikes for initialization
rez                = kilosort('fitTemplates', rez, DATA, uproj);  % fit templates iteratively
rez                = kilosort('fullMPMU', rez, DATA);% extract final spike times (overlapping extraction)

if fMerge_post
    rez = kilosort('merge_posthoc2', rez); %ask whether to merge or not
end

% save python results file for Phy
if fSavePhy
    try
        kilosort('rezToPhy', rez, fpath); %path to npy2mat needed
    catch
        disperr_();
    end
end

runtime_ksort = toc(runtime_ksort);
fprintf('\tkilosort took %0.1fs for %s\n', runtime_ksort, P.vcFile_prm);

% output kilosort result
S_ksort = struct('rez', rez, 'P', P, 'runtime_ksort', runtime_ksort);
struct_save_(S_ksort, strrep(vcFile_prm, '.prm', '_ksort.mat'), 1);
end %func


%--------------------------------------------------------------------------
function struct_save_(S, vcFile, fVerbose)
% 7/13/17 JJJ: Version check routine
if nargin<3, fVerbose = 0; end
if fVerbose
    fprintf('Saving a struct to %s...\n', vcFile); t1=tic;
end
version_year = version('-release');
version_year = str2double(version_year(1:end-1));
if version_year >= 2017
    save(vcFile, '-struct', 'S', '-v7.3', '-nocompression'); %faster    
else
%     disp('Saving with -nocompression flag failed. Trying without compression');
    save(vcFile, '-struct', 'S', '-v7.3');
end
if fVerbose
    fprintf('\ttook %0.1fs.\n', toc(t1));
end
end %func


%--------------------------------------------------------------------------
function export_imec_sync_(vcFiles_prm)
% sync output for IMEC (export the last channel entries)
try
    P = file2struct_(vcFiles_prm);
    vcFile_bin = P.vcFile;    
    fid = fopen(vcFile_bin, 'r');
    vnSync = fread(fid, inf, '*uint16'); fclose(fid);
    vnSync = vnSync(P.nChans:P.nChans:end); %subsample sync channel
    assignWorkspace_(vnSync);
catch
    fprintf(2, 'error exporting sync: %s\n', lasterr());
end
end


%--------------------------------------------------------------------------
function plot_activity_(P) % single column only
% plot activity as a function of depth and time
global fDebug_ui

tbin = 10; %activity every 10 sec
% plot activity as a function of time
% vcFile_evt = subsFileExt(P.vcFile_prm, '_evt.mat');
S0 = load_cached_(P, 0); %do not load waveforms
nSites = numel(P.viSite2Chan);
% tdur = max(cell2mat_(cellfun(@(x)double(max(x)), Sevt.cviSpk_site, 'UniformOutput', 0))) / P.sRateHz;
tdur = double(max(S0.viTime_spk)) / P.sRateHz; % in sec
nTime = ceil(tdur / tbin);

mrAmp90 = zeros(nTime, nSites);
lim0 = [1, tbin * P.sRateHz];
for iSite=1:nSites
    viSpk1 = find(S0.viSite_spk == iSite);
    vrAmp_spk1 = S0.vrAmp_spk(viSpk1); % % S0.cvrSpk_site{iSite};  %spike amplitude   
    if isempty(vrAmp_spk1), continue; end
    viTime_spk1 = S0.viTime_spk(viSpk1);
    for iTime=1:nTime
        lim1 = lim0 + (iTime-1) * lim0(2);
        vrAmp_spk11 = vrAmp_spk1(viTime_spk1 >= lim1(1) & viTime_spk1 <= lim1(2));
        if isempty(vrAmp_spk11),  continue; end
        mrAmp90(iTime, iSite) = quantile(abs(vrAmp_spk11), .9);
    end
end %for
mrAmp90=mrAmp90';

vlSite_left = P.mrSiteXY(:,1) == 0;
vrSiteY = P.mrSiteXY(:,2);
hFig = create_figure_('FigActivity', [0 0 .5 1], P.vcFile_prm, 1, 1);
subplot 121; imagesc(mrAmp90(vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites'); title('Left edge sites');
subplot 122; imagesc(mrAmp90(~vlSite_left, :), 'XData', (1:nTime) * tbin, 'YData', vrSiteY(~vlSite_left)); axis xy; xlabel('Time'); ylabel('Sites'); title('Right edge sites');

[~, iSite_center] = max(mean(mrAmp90,2)); 
viSiteA = iSite_center + [-2:2]; %look neighbors
mrAmp90a = mrAmp90(viSiteA, :);
vrCentroid = bsxfun(@rdivide, sum(bsxfun(@times, mrAmp90a.^2, vrSiteY(viSiteA))), sum(mrAmp90a.^2));
hold on; plot((1:nTime) * tbin, vrCentroid, 'r');

% if get_set_([], 'fDebug_ui', 0), close(hFig); end
if fDebug_ui==1, close(hFig); end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: multiple targeting copy file. Tested
function copyfile_(csFiles, vcDir_dest)
% copyfile_(vcFile, vcDir_dest)
% copyfile_(csFiles, vcDir_dest)
% copyfile_(csFiles, csDir_dest)

% Recursion if cell is used
if iscell(vcDir_dest)
    csDir_dest = vcDir_dest;
    for iDir = 1:numel(csDir_dest)
        try
            copyfile_(csFiles, csDir_dest{iDir});
        catch
            disperr_();
        end
    end
    return;
end

if ischar(csFiles), csFiles = {csFiles}; end
for iFile=1:numel(csFiles)
    vcPath_from_ = csFiles{iFile};
    if exist(vcPath_from_, 'dir') == 7        
        [vcPath_,~,~] = fileparts(vcPath_from_);
        vcPath_from_ =  sprintf('%s%s*', vcPath_, filesep());
        vcPath_to_ = sprintf('%s%s%s%s', vcDir_dest, filesep(), vcPath_, filesep());
        if exist(vcPath_to_, 'dir') ~= 7, mkdir(vcPath_to_); end    
        disp([vcPath_from_, '; ', vcPath_to_]);
    else
        vcPath_to_ = vcDir_dest;
        if exist(vcPath_to_, 'dir') ~= 7
            mkdir(vcPath_to_); 
            disp(['Created a folder ', vcPath_to_]);
        end
    end
    try   
        vcEval1 = sprintf('copyfile ''%s'' ''%s'' f;', vcPath_from_, vcPath_to_);
        eval(vcEval1);
        fprintf('\tCopied %s to %s\n', vcPath_from_, vcPath_to_);
    catch
        fprintf(2, '\tFailed to copy %s\n', vcPath_from_);
    end
end
end %func


%--------------------------------------------------------------------------
% 9/17/17 JJJ: nGroups fixed, can have less than 7 clusters
function vhPlot = plot_group_(hAx, mrX, mrY, varargin)
nGroups = min(7, size(mrX,2));
vhPlot = zeros(nGroups, 1);
hold(hAx, 'on');
for iGroup=1:nGroups
    vrX1 = mrX(:, iGroup:nGroups:end);
    vrY1 = mrY(:, iGroup:nGroups:end);
    vhPlot(iGroup) = plot(hAx, vrX1(:), vrY1(:), varargin{:});
end
end %func


%--------------------------------------------------------------------------
function plot_update_(vhPlot, mrX, mrY)
nPlots = numel(vhPlot);
for iPlot=1:nPlots
    vrX1 = mrX(:, iPlot:nPlots:end);
    vrY1 = mrY(:, iPlot:nPlots:end);
    set(vhPlot(iPlot), 'XData', vrX1(:), 'YData', vrY1(:));
end
end %func


%--------------------------------------------------------------------------
function varargout = gather_(varargin)
for i=1:nargin
    try
        varargout{i} = gather(varargin{i});
    catch
        varargout{i} = [];
    end
end
end %func


%--------------------------------------------------------------------------
function vc = set_(vc, varargin)
% Set handle to certain values
% set_(S, name1, val1, name2, val2)

if isempty(vc), return; end
if isstruct(vc)
    for i=1:2:numel(varargin)        
        vc.(varargin{i}) = varargin{i+1};
    end
    return;
end
if iscell(vc)
    for i=1:numel(vc)
        try
            set(vc{i}, varargin{:});
        catch
        end
    end
elseif numel(vc)>1
    for i=1:numel(vc)
        try
            set(vc(i), varargin{:});
        catch
        end
    end
else
    try
        set(vc, varargin{:});
    catch
    end 
end
end %func


%--------------------------------------------------------------------------
function vr = tnWav2uV_(vn, P)
% Integrate the waveform

if nargin<2, P = get0_('P'); end
vr = bit2uV_(vn, P);
% if P.nDiff_filt>0, 
vr = meanSubt_(vr); 
% end
end %func


%--------------------------------------------------------------------------
function export_spkwav_(P, vcArg2, fDiff)
% Export spike waveforms organized by clusters
% export_spkwav_(P)
% export_spkwav_(P, viClu)
global fDebug_ui

if nargin<2, vcArg2 = ''; end
if nargin<3, fDiff = 0; end
if isempty(vcArg2)
    viClu = [];
    fPlot = 0;
else
    viClu = str2num(vcArg2);
    if isempty(viClu), fprintf(2, 'Invalid Cluster #: %s', vcArg2); end
    fPlot = numel(viClu)==1;
end

% Load data
S0 = load_cached_(P);
if isempty(S0), fprintf(2, 'Not clustered yet'); return; end
if ~isfield(S0, 'S_clu'), fprintf(2, 'Not clustered yet'); return; end
S_clu = S0.S_clu;
% global tnWav_spk

% Collect waveforms by clusters
ctrWav_clu = cell(1, S_clu.nClu);
miSite_clu = P.miSites(:,S_clu.viSite_clu);
fprintf('Collecting spikes from clusters\n\t'); t1=tic;
if isempty(viClu), viClu = 1:S_clu.nClu; end
for iClu = viClu
    tnWav_clu1 = tnWav_spk_sites_(S_clu.cviSpk_clu{iClu}, miSite_clu(:,iClu));
    if fDiff
        ctrWav_clu{iClu} = tnWav_clu1;
    else
        ctrWav_clu{iClu} = tnWav2uV_(tnWav_clu1);
    end
%     ctrWav_clu{iClu} = (cumsum(bit2uV_(meanSubt_(tnWav_clu1))));
%     ctrWav_clu{iClu} = (bit2uV_(meanSubt_(tnWav_clu1)));
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t1));

if fPlot
    iClu = viClu;
    nT_spk = (diff(P.spkLim)+1);
    nSpk1 = S_clu.vnSpk_clu(iClu);
    hFig = create_figure_(sprintf('Fig_clu%d', iClu), [0 0 .5 1], P.vcFile_prm, 1, 1);
    multiplot([], P.maxAmp, [], ctrWav_clu{iClu}, miSite_clu(:,iClu));
    xlabel('Spike #'); ylabel('Site #'); 
    set(gca, 'YTick', range_(miSite_clu(:,iClu)), ...
        'XTick', (1:nSpk1) * nT_spk - P.spkLim(2), 'XTickLabel', 1:nSpk1); 
    axis([0, nSpk1*nT_spk+1, (limit_(miSite_clu(:,iClu)) + [-1,1])]);
    grid on;
    title(sprintf('Cluster%d (n=%d), Scale=%0.1f uV', iClu, nSpk1, P.maxAmp)); 
    mouse_figure(hFig);
    eval(sprintf('trWav_clu%d = ctrWav_clu{iClu};', iClu));
    eval(sprintf('viSites_clu%d = miSite_clu(:,iClu);', iClu));
    eval(sprintf('assignWorkspace_(trWav_clu%d, viSites_clu%d);', iClu, iClu));
%     if get_set_([], 'fDebug_ui', 0), close_(hFig); end
    if fDebug_ui==1, close_(hFig); end
else
    assignWorkspace_(ctrWav_clu, miSite_clu);
end
end %func


%--------------------------------------------------------------------------
function range1 = range_(vn)
vn = vn(:);
range1 = min(vn):max(vn);
end


%--------------------------------------------------------------------------
function limit1 = limit_(vn)
vn = vn(:);
limit1 = [min(vn), max(vn)];
end


%--------------------------------------------------------------------------
function try_eval_(vcEval1, fVerbose)
if nargin<2, fVerbose = 1; end
try 
    eval(vcEval1); 
    fprintf('\t%s\n', vcEval1);  
catch
    if fVerbose
        fprintf(2, '\tError evaluating ''%s''\n', vcEval1);  
    end
end
end %func


%--------------------------------------------------------------------------
function export_spkamp_(P, vcArg2)
% export spike amplitudes (Vpp, Vmin, Vmax) in uV to workspace, organize by clusters

if nargin<2, vcArg2 = ''; end
if isempty(vcArg2)
    viClu = [];
    fSingleUnit = 0;
else
    viClu = str2num(vcArg2);
    if isempty(viClu), fprintf(2, 'Invalid Cluster #: %s', vcArg2); end
    fSingleUnit = numel(viClu)==1;
end

% Load data
S0 = load_cached_(P);
if isempty(S0), fprintf(2, 'Not clustered yet'); return; end
if ~isfield(S0, 'S_clu'), fprintf(2, 'Not clustered yet'); return; end
S_clu = S0.S_clu;
% global tnWav_spk

% Collect waveforms by clusters
[cmrVpp_clu, cmrVmin_clu, cmrVmax_clu] = deal(cell(1, S_clu.nClu));
miSite_clu = P.miSites(:,S_clu.viSite_clu);
fprintf('Calculating spike amplitudes from clusters\n\t'); t1=tic;
if isempty(viClu), viClu = 1:S_clu.nClu; end
for iClu = viClu
    trWav_clu1 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu}, miSite_clu(:,iClu)));
    cmrVmin_clu{iClu} = shiftdim(min(trWav_clu1));
    cmrVmax_clu{iClu} = shiftdim(max(trWav_clu1));
    cmrVpp_clu{iClu} = cmrVmax_clu{iClu} - cmrVmin_clu{iClu};
    fprintf('.');
end
fprintf('\n\ttook %0.1fs\n', toc(t1));

% export to workspace
if fSingleUnit
    iClu = viClu;
    eval(sprintf('mrVmax_clu%d = cmrVmax_clu{iClu};', iClu));
    eval(sprintf('mrVmin_clu%d = cmrVmin_clu{iClu};', iClu));
    eval(sprintf('mrVpp_clu%d = cmrVpp_clu{iClu};', iClu));
    eval(sprintf('viSites_clu%d = miSite_clu(:,iClu);', iClu));
    eval(sprintf('assignWorkspace_(mrVmax_clu%d, mrVmin_clu%d, mrVpp_clu%d, viSites_clu%d);', iClu, iClu, iClu, iClu)); 
else
    assignWorkspace_(cmrVpp_clu, cmrVmax_clu, cmrVmin_clu, miSite_clu);
end
end %func


%--------------------------------------------------------------------------
function csFiles_bin = dir_files_(csFile_merge, vcFile_txt, file_sort_merge)
% Display binary files
if nargin<2, vcFile_txt=''; end
if nargin<3, file_sort_merge = []; end
if isempty(file_sort_merge)
    file_sort_merge = 1; 
elseif ischar(file_sort_merge)
    file_sort_merge = str2double(file_sort_merge);
end

if isempty(csFile_merge)
    fprintf('No files to merge ("csFile_merge" is empty).\n');
    csFiles_bin = {}; return; 
end
% fprintf('Listing files to merge ("csFile_merge"):\n');
csFiles_bin = filter_files_(csFile_merge, file_sort_merge);
if nargout==0
%     arrayfun(@(i)fprintf('%d: %s\n', i, csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
    arrayfun(@(i)fprintf('%s\n', csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
end
if ~isempty(vcFile_txt)
    cellstr2file_(vcFile_txt, csFiles_bin);
    fprintf('%s is created.\n', vcFile_txt);
    edit(vcFile_txt);
end
end %func


%--------------------------------------------------------------------------
function commit_jrc3_(S_cfg, vcDir, fZipFile)
global fDebug_ui

if nargin<1, S_cfg = read_cfg_(); end 
if nargin<2, vcDir = []; end
if nargin<3, fZipFile = 1; end

fDebug_ui = 0;
set0_(fDebug_ui);
if isempty(vcDir)
    vcDir = [S_cfg.path_dropbox3, filesep()];
end

disp(['Commiting to ', vcDir]);
for i=1:2
    disp(['Deleting files in ' vcDir]);
    delete_([vcDir, '*']);
%     delete_([vcDir, '\kilosort\*']);
    pause(.5);
end

copyfile_(S_cfg.sync_list_ver3, vcDir); %destination root
copyfile_(S_cfg.csFiles_cu3, vcDir); %destination root
copyfile_(S_cfg.csFiles_ptx3, vcDir); %destination root
% mkdir_([vcDir, 'kilosort']);
% copyfile_('./kilosort/*', [vcDir, 'kilosort']); %destination root
% copyfile_(S_cfg.csFiles_sample, vcDir);

% Zip 
% delete_([vcDir, 'jrc3.zip']); %delete previously zipped file
if fZipFile
    hMsg = msgbox_(sprintf('Archiving to %s', [vcDir, 'jrc3.zip']));
    t1 = tic;
    [csFiles_jrc3_full, csFiles_jrc3] = dir_([vcDir, '*'], 'jrc3.zip');
    zip([vcDir, 'jrc3.zip'], csFiles_jrc3, vcDir);
    fprintf('Zip file creation took %0.1f\n', toc(t1));
    close_(hMsg);    
    msgbox_('Update the Dropbox link for www.jrclust.org');
end
end %func


%--------------------------------------------------------------------------
function mkdir_(vcDir)
% make only if it doesn't exist. provide full path for dir
if ~exist(vcDir, 'dir'), mkdir(vcDir); end
end %func


%--------------------------------------------------------------------------
function [csFiles_full, csFiles] = dir_(vcFilter_dir, csExcl)
% return name of files full path, exclude files
if nargin>=2
    if ischar(csExcl), csExcl = {csExcl}; end
    csExcl = union(csExcl, {'.', '..'}); 
else
    csExcl = [];
end
csFiles = dir(vcFilter_dir);
csFiles  = {csFiles.('name')};
csFiles = setdiff(csFiles, csExcl);
[vcDir, ~, ~] = fileparts(vcFilter_dir);
if isempty(vcDir), vcDir='.'; end
csFiles_full = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
function delete_(csFiles)
if ischar(csFiles), csFiles = {csFiles}; end
for i=1:numel(csFiles)
    try
        if iscell(csFiles)
            delete(csFiles{i});
        else
            delete(csFiles(i));
        end
    catch
%         disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
function export_fet_(P)
% export feature matrix to workspace
S0 = load(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
trFet = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
%mrFet2 = load_bin_(strrep(P.vcFile_prm, '.prm', '_fet2.bin'), 'single', S0.dimm_fet);
%miFet_sites = load_bin_(strrep(P.vcFile_prm, '.prm', '_fet_sites.bin'), 'int32', S0.dimm_fet_sites);
%assignWorkspace_(mrFet, miFet_sites);
assignWorkspace_(trFet);
end %func


%--------------------------------------------------------------------------
function varargout = getfield_(S, varargin)
for iField = 1:numel(varargin)
    if isfield(S, varargin{iField})
        varargout{iField} = getfield(S, varargin{iField});
    else
        varargout{iField} = [];
    end
end
end %func


%--------------------------------------------------------------------------
function [mrFet, vrY_pre, vrY_post] = drift_correct_(mrFet, P) % drift correction
% Manually correct drift using {[t1,t2,d1,d2], ...} format in cvrDepth_drift
% mrFet last column contains the y coordinate
fPlot_drift = 1;
if isempty(getfield_(P, 'cvrDepth_drift')), return; end
S0 = get0_();
fprintf('Correcting drift\n\t'); t1=tic;
vrY_pre = mrFet(end,:);
vrY_post = vrY_pre;
for iDrift = 1:numel(P.cvrDepth_drift)
    tlim_dlim1 = P.cvrDepth_drift{iDrift};
    tlim1 = tlim_dlim1(1:2) * P.sRateHz;
    dlim1 = tlim_dlim1(3:4);
    viSpk1 = find(S0.viTime_spk >= tlim1(1) & S0.viTime_spk < tlim1(2));
    if isempty(viSpk1), continue; end
    viTime_spk1 = S0.viTime_spk(viSpk1);
    if diff(dlim1) == 0
        vrY_post(viSpk1) = vrY_pre(viSpk1) + dlim1(1); % use single depth correction factor 
    else        
        vrY_post(viSpk1) = vrY_pre(viSpk1) + interp1(tlim1, dlim1, viTime_spk1, 'linear', 'extrap');
    end    
    fprintf('.');
end 
mrFet(end,:) = vrY_post;
fprintf('\n\tDrift correction took %0.1fs\n', toc(t1));

if fPlot_drift
    [viTime_spk, vrAmp_spk] = get0_('viTime_spk', 'vrAmp_spk');
    viSpk = find(vrAmp_spk < median(vrAmp_spk)); %pick more negative
%     viSpk1 = viSpk1(1:2:end); %plot every other
    figure; hold on; set(gcf,'Color','w');
    ax(1)=subplot(121); 
    plot(viTime_spk(viSpk), vrY_pre(viSpk), 'go', 'MarkerSize', 2); title('before correction'); grid on;
    ax(2)=subplot(122); 
    plot(viTime_spk(viSpk), vrY_post(viSpk), 'bo', 'MarkerSize', 2); title('after correction'); grid on;
    linkaxes(ax,'xy');
end
end %func


%--------------------------------------------------------------------------
function mrDist_clu = S_clu_cov_(S_clu, P)
% merge clusters based on the spike waveforms
% nPc = 3;
% trCov = tr2cov_(S_clu.trWav_spk_clu);
% trWav_clu = meanSubt_(S_clu.trWav_spk_clu);
trWav_clu = meanSubt_(S_clu.trWav_raw_clu); %uses unfiltered waveform
% pairwise similarity computation
% maxSite = P.maxSite_merge;
maxSite = ceil(P.maxSite/2);
% maxSite = P.maxSite;
mrDist_clu = nan(S_clu.nClu);
vrVar_clu = zeros(S_clu.nClu, 1, 'single');
% miSites = P.miSites(1:P.maxSite_merge*2+1, :);
for iClu2 = 1:S_clu.nClu        
    iSite2 = S_clu.viSite_clu(iClu2);
%     viSite2 = miSites(:,iSite2);
    viClu1 = find(abs(S_clu.viSite_clu - iSite2) <= maxSite);
    viClu1(viClu1 == iClu2) = [];
    viClu1 = viClu1(:)';
%     mrCov2 = eigvec_(trCov(:,:,iClu2), nPc);    
    mrWav_clu2 = trWav_clu(:,:,iClu2);
    [~,~,var2] = pca(mrWav_clu2); var2 = var2(1)/sum(var2);
    vrVar_clu(iClu2) = var2;
    for iClu1 = viClu1           
%         mrCov1 = eigvec_(trCov(:,:,iClu1), nPc);
        mrWav_clu1 = trWav_clu(:,:,iClu1);
%         mrDist_clu(iClu1, iClu2) = mean(abs(mean(mrCov1 .* mrCov2))); %cov_eig_(mrCov1, mrCov2, nPc);
        mrDist_clu(iClu1, iClu2) = cluWav_dist_(mrWav_clu1, mrWav_clu2);
    end
end
end %func


%--------------------------------------------------------------------------
function trCov = tr2cov_(trWav)
[nT, nSites, nClu] = size(trWav);
trCov = zeros(nT, nT, nClu, 'like', trWav);
trWav = meanSubt_(trWav);
for iClu=1:nClu
    mrCov1 = trWav(:,:,iClu); %mean waveform covariance
    mrCov1 = mrCov1 * mrCov1';    
    trCov(:,:,iClu) = mrCov1;
end %for
end %func


%--------------------------------------------------------------------------
function dist12 = cov_eig_(mrCov1, mrCov2, nPc)
d1 = cov2var_(mrCov1, nPc); 
d2 = cov2var_(mrCov2, nPc);
d12 = cov2var_(mrCov1 + mrCov2, nPc);
dist12 = d12 / ((d1+d2)/2);
end %func


%--------------------------------------------------------------------------
function [vrD1, mrPv1] = cov2var_(mrCov1, nPc)
[mrPv1, vrD1] = eig(mrCov1); 
vrD1 = cumsum(flipud(diag(vrD1))); 
vrD1 = vrD1 / vrD1(end);
if nargin>=2
    vrD1 = vrD1(nPc);
end
end


%--------------------------------------------------------------------------
function mrCov2 = eigvec_(mrCov2, nPc)
[mrCov2,~] = eig(mrCov2); 
mrCov2 = mrCov2(:, end-nPc+1:end);
mrCov2 = bsxfun(@minus, mrCov2, mean(mrCov2));
mrCov2 = bsxfun(@rdivide, mrCov2, std(mrCov2));
end %func


%--------------------------------------------------------------------------
function corr12 = cluWav_dist_(mrWav_clu1, mrWav_clu2)
nPc = 3;
viDelay = 0;

[~, mrPv1, vrL1] = pca(mrWav_clu1, 'NumComponents', nPc);
[~, mrPv2, vrL2] = pca(mrWav_clu2, 'NumComponents', nPc);
% viDelay = -4:4; %account for time delay

if viDelay==0
    corr12 = mean(abs(diag(corr_(mrPv1, mrPv2))));
else
    vrDist12 = zeros(size(viDelay), 'single');
    for iDelay=1:numel(viDelay)
        mrPv2a = shift_mr_(mrPv2, viDelay(iDelay));
        vrDist12(iDelay) = mean(abs(diag(corr_(mrPv1, mrPv2a))));
    end
    corr12 = max(vrDist12);
end
% dist12 = abs(corr_(mrPv1, mrPv2));
% [~,~,vrL12] = pca([mrWav_clu1, mrWav_clu2]);
% nPc = 1;
% dist12 = sum(vrL12(1:nPc)) / sum(vrL12);
% 
% [~,~,vrL1] = pca(mrWav_clu1);
% dist1 = sum(vrL1(1:nPc)) / sum(vrL1);
% % dist1 = vrL(1) / sum(vrL);
% % 
% [~,~,vrL2] = pca(mrWav_clu2);
% dist2 = sum(vrL2(1:nPc)) / sum(vrL2);
% % dist2 = vrL(1) / sum(vrL);
% % 
% dist12 = dist12 / ((dist1 + dist2)/2);

% mrCov12 = [mrWav_clu1, mrWav_clu2]; 
% mrCov12 = mrCov12 * mrCov12';
% [~,vrD] = eig([mrWav_clu1, mrWav_clu2]);
% dist12 = vrD(end) / sum(vrD);

% nPc = 3;
% mrPc12 = eigvec_(mrCov12 * mrCov12', nPc);
% mrPc1 = eigvec_(mrWav_clu1 * mrWav_clu1', nPc);
% mrPc2 = eigvec_(mrWav_clu2 * mrWav_clu2', nPc);
% dist12 = mean(abs(mean(mrPc2 .* mrPc1)));
end %func


%--------------------------------------------------------------------------
function mrDist_clu = S_clu_spkcov_(S_clu, P)
% compute covariance from each spikes belonging to clusters
global tnWav_raw %tnWav_spk 
MAX_SAMPLE = 1000;
nPc = 3;
viSite_spk = get0_('viSite_spk');
maxSite = ceil(P.maxSite/2);
mrDist_clu = nan(S_clu.nClu);

% compute a cell of mrPv
[cmrPv_spk_clu, cmrPv_raw_clu] = deal(cell(S_clu.nClu, 1));
for iClu = 1:S_clu.nClu
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    viSpk_clu1 = viSpk_clu1(viSite_spk(viSpk_clu1) == S_clu.viSite_clu(iClu));
    viSpk_clu1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
%     cmrPv_spk_clu{iClu} = tn2pca_spk_(tnWav_spk(:,:,viSpk_clu1), nPc);
    cmrPv_raw_clu{iClu} = tn2pca_spk_(tnWav_raw(:,:,viSpk_clu1), nPc);
end %for

for iClu2 = 1:S_clu.nClu        
    iSite2 = S_clu.viSite_clu(iClu2);
    viClu1 = find(abs(S_clu.viSite_clu - iSite2) <= maxSite);
    viClu1(viClu1 == iClu2) = [];
    viClu1 = viClu1(:)';
%     mrPv_clu2 = cmrPv_clu_raw{iClu2};
    for iClu1 = viClu1           
        vrCorr_raw = diag(corr_(cmrPv_raw_clu{iClu2}, cmrPv_raw_clu{iClu1}));
%         vrCorr_spk = diag(corr_(cmrPv_spk_clu{iClu2}, cmrPv_spk_clu{iClu1}));
%         mrDist_clu(iClu1, iClu2) = mean(abs([vrCorr_raw; vrCorr_spk]));
%         mrDist_clu(iClu1, iClu2) = mean(abs([vrCorr_raw; vrCorr_spk]));
%         mrPv_clu1 = cmrPv_clu_raw{iClu1};        
        mrDist_clu(iClu1, iClu2) = mean(abs(vrCorr_raw));
    end
end %for
end %func


%--------------------------------------------------------------------------
function mrPv = tn2pca_cov_(tn, nPc)
% calc cov by averaging and compute pca
[nT, nC, nSpk] = size(tn);
mrCov1 = zeros(nT, nT, 'double');
tr = meanSubt_(single(tn));
for iSpk = 1:nSpk
    mr1 = tr(:,:,iSpk);
    mrCov1 = mrCov1 + double(mr1 * mr1');
end
[mrPv, vrL] = eig(mrCov1);
mrPv = mrPv(:, end:-1:end-nPc+1);
end %func


%--------------------------------------------------------------------------
function mrPv = tn2pca_spk_(tn, nPc)
% calc cov by averaging and compute pca
[nT, nC, nSpk] = size(tn);
% mrCov1 = zeros(nT, nT, 'single');
mr = single(reshape(tn, nT,[]));
mr = meanSubt_(mr);
[~, mrPv] = pca(mr, 'NumComponents', nPc);
% tr = meanSubt_(single(tn));
% for iSpk = 1:nSpk
%     mr1 = tr(:,:,iSpk);
%     mrCov1 = mrCov1 + mr1 * mr1';
% end
% [mrPv, vrL] = eig(mrCov1);
% mrPv = mrPv(:, end:-1:end-nPc+1);
end %func


%--------------------------------------------------------------------------
function mr = shift_vr_(vr, vn_shift)
n = numel(vr);
% mr = zeros(n, numel(vn_shift), 'like', vr);
mi = bsxfun(@plus, (1:n)', vn_shift(:)');
mi(mi<1) = 1;
mi(mi>n) = n;
mr = vr(mi);
end %func


%--------------------------------------------------------------------------
function [iClu, hPlot] = plot_tmrWav_clu_(S0, iClu, hPlot, vrColor)
[S_clu, P] = getfield_(S0, 'S_clu', 'P');
[hFig, S_fig] = get_fig_cache_('FigWav');
% figure(hFig);
% if isempty(get_(hPlot, 'Parent')), delete_(hPlot); hPlot = []; end
% if ~isvalid(hPlot), hPlot = []; end
if isempty(hPlot)
    hPlot = plot(nan, nan, 'Color', vrColor, 'LineWidth', 2, 'Parent', S_fig.hAx);
end
if P.fWav_raw_show
    mrWav_clu1 = S_clu.tmrWav_raw_clu(:,:, iClu);
else
    mrWav_clu1 = S_clu.tmrWav_clu(:,:, iClu);
end
multiplot(hPlot, S_fig.maxAmp, wav_clu_x_(iClu, P), mrWav_clu1);
uistack_(hPlot, 'top');
end %func


%--------------------------------------------------------------------------
function mrWav = mean_tnWav_raw_(tnWav, P)
mrWav = meanSubt_(mean(single(tnWav),3)) * P.uV_per_bit;
end %func


%--------------------------------------------------------------------------
function tr = raw2uV_(tnWav_raw, P)
tr = meanSubt_(single(tnWav_raw) * P.uV_per_bit);
end %func


%--------------------------------------------------------------------------
function raw_waveform_(hMenu)
global tnWav_spk tnWav_raw

[P, S_clu] = get0_('P', 'S_clu');
if get_(P, 'fWav_raw_show')
    P.fWav_raw_show = 0;
else
    P.fWav_raw_show = 1;
end

if isempty(get_(S_clu, 'tmrWav_raw_clu'))    
    % check for backward compatibliity
    viSite_spk = get0_('viSite_spk');
    S_clu = clu2wav_(S_clu, viSite_spk, tnWav_spk, tnWav_raw);
    S0 = set0_(P, S_clu);
else
    S0 = set0_(P);
end
set(hMenu, 'Checked', ifeq_(P.fWav_raw_show, 'on', 'off'));
% redraw windows
plot_FigWav_(S0);
button_CluWav_simulate_(S0.iCluCopy, S0.iCluPaste, S0);
end %func


%--------------------------------------------------------------------------
function [mrPv1, mrPv2] = pca_pv_spk_(viSpk1, viSites1, tnWav_spk1)
% if viSite not found just set to zero
% error('not implemented yet');
% global mrPv_global
if nargin<3
    tnWav_spk1 = permute(tnWav_spk_sites_(viSpk1, viSites1), [1,3,2]);
end
nT = size(tnWav_spk1, 1);
nSites = numel(viSites1);
mrPv_global = get0_('mrPv_global');
if ~isempty(mrPv_global)
    % show global pca
    mrPv1 = repmat(mrPv_global(:,1), [1, nSites]);
    mrPv2 = repmat(mrPv_global(:,2), [1, nSites]);
else
    % show site pca
    [mrPv1, mrPv2] = deal(zeros(nT, nSites, 'single'));
    for iSite1=1:nSites
        [mrPv1(:,iSite1), mrPv2(:,iSite1)] = pca_pv_(tnWav_spk1(:,:,iSite1));
    end %for
end
end %func


%--------------------------------------------------------------------------
function [vrPv1, vrPv2] = pca_pv_(mr1)
MAX_SAMPLE = 10000; %for pca
mrCov = subsample_mr_(mr1, MAX_SAMPLE, 2);
mrCov = meanSubt_(single(mrCov));
mrCov = mrCov * mrCov';
[mrPv1, vrD1] = eig(mrCov); 
vrPv1 = mrPv1(:,end);
vrPv2 = mrPv1(:,end-1);
% iMid = 1-P.spkLim(1);
% vrPv1 = ifeq_(vrPv1(iMid)>0, -vrPv1, vrPv1);
% vrPv2 = ifeq_(vrPv2(iMid)>0, -vrPv2, vrPv2);
end %func


%--------------------------------------------------------------------------
function [mrPc1, mrPc2, mrPv1, mrPv2] = pca_pc_spk_(viSpk1, viSites1, mrPv1, mrPv2)
% varargin: mrPv
% varargout: mrPc
% project
%[mrPc1, mrPc2, mrPv1, mrPv2] = pca_pc_spk_(viSpk1, viSites1)
%[mrPc1, mrPc2] = pca_pc_spk_(viSpk1, viSites1, mrPv1, mrPv2)
nSites1 = numel(viSites1);
tnWav_spk1 = permute(tnWav_spk_sites_(viSpk1, viSites1), [1,3,2]);
if nargin<3
    [mrPv1, mrPv2] = pca_pv_spk_(viSpk1, viSites1, tnWav_spk1);
end

dimm1 = size(tnWav_spk1); %nT x nSpk x nChan
[mrPc1, mrPc2] = deal(zeros(dimm1(2), nSites1, 'single'));
try
    for iSite1=1:nSites1
        mrWav_spk1 = meanSubt_(single(tnWav_spk1(:,:,iSite1)));
        mrPc1(:,iSite1) = (mrPv1(:,iSite1)' * mrWav_spk1)';
        mrPc2(:,iSite1) = (mrPv2(:,iSite1)' * mrWav_spk1)';
    end %for
catch
    disperr_();
end
mrPc1 = (mrPc1') / dimm1(1);
mrPc2 = (mrPc2') / dimm1(1);
end %func


%--------------------------------------------------------------------------
function vc = if_on_off_(vc, cs)
if ischar(cs), cs = {cs}; end
vc = ifeq_(ismember(vc, cs), 'on', 'off');
end %func


%--------------------------------------------------------------------------
function proj_view_(hMenu)
P = get0_('P');
vcFet_show = lower(get(hMenu, 'Label'));
switch vcFet_show
    case {'vpp', 'pca', 'cov'}, P.vcFet_show = vcFet_show;
    otherwise, error('proj_view_:Invalid vcFet_show');
end
vhMenu = hMenu.Parent.Children;
for iMenu=1:numel(vhMenu)
    vhMenu(iMenu).Checked = if_on_off_(vhMenu(iMenu).Label, vcFet_show);
end
% auto-scale the view
set0_(P);
auto_scale_proj_time_();
end %func


%--------------------------------------------------------------------------
function P = struct_default_(P, csName, def_val)
% Set to def_val if empty or field does not exist
% set the field(s) to default val

if ischar(csName), csName = {csName}; end
for iField = 1:numel(csName)
    vcName = csName{iField};
    if ~isfield(P, vcName)
        P.(vcName) = def_val;
    elseif isempty(P.(vcName))
        P.(vcName) = def_val;
    end
end
end %func


% --------------------------------------------------------------------------
% function [S_clu, nMerges_clu] = S_clu_isosplit_(S_clu)
% % Jeremy Maglund inspired
% error('not implemented yet');
% % simple merge only, no cluster refinement
% % average cluster waveform based merging only, no individual spike split @TODO
% 
% % global tnWav_spk
% 
% nRepeat = 3;
% 
% P = get0_('P');
% % Build a distance matrix and identify the clostest pairs
% % comparisons_made = ;
% % final_pass = false;
% % labels = S_clu.viClu;
% % Kmax = S_clu.nClu;
% mrFet_clu = compute_centers_(S_clu); %nDim x nClu
% mrDist_clu = make_dists_matrix_(mrFet_clu, P); % Kmax x Kmax, inf diagonal
% 
% % preserve cluster index while merging and compact it later
% for iRepeat = 1:nRepeat
% %     active_labels = unique(labels);
% %     mrFet_mean_clu = fet_mean_clu_(S_clu);
% %     mrDist_clu = dist_clu_(S_clu, mrFet_mean_clu); % 
%     [inds1, inds2] = get_pairs_to_compare_(mrDist_clu);
%     if isempty(inds1), break; end
% 
%     [S_clu, mrFet_clu, vi_inds_merged] = compare_pairs_(S_clu, mrFet_clu, inds1, inds2);    
%     if isempty(vi_inds_merged), return; end %@TODO: may have to repartition later
%     
%     % update changed clusters and distance compared    
%     mrDist_clu = update_dists_matrix_(mrDist_clu, mrFet_clu, inds1, inds2, vi_inds_merged);    
% end %for
% 
% % [miClu_pairs, mrDist_clu] = clu_dist();
% % S_iso = struct('isocut_threshold', 1, 'refine_clusters', 1, 'max_iterations_per_pass', 500, 'whiten_cluster_pairs', 1, 'prevent_merge', false);
% % 
% % mrFet_center_clu = compute_centers_(mrFet, S_clu.viClu); 
% % final_pass = false;
% % comparisons_made=zeros(Kmax,Kmax);
% % while 1
% %     something_merged = false;
% %     clusters_changed_vec_in_pass=zeros(1,Kmax);
% %     for iLoop = 1:S_iso.max_iterations_per_pass
% %         
% %     end %for
% %     % zero out the comparisons made matrix only for those that have changed
% %     clusters_changed=find(clusters_changed_vec_in_pass);
% %     for j=1:length(clusters_changed)
% %         comparisons_made(clusters_changed(j),:)=0;
% %         comparisons_made(:,clusters_changed(j))=0;
% %     end;
% %     
% %     if (something_merged), final_pass=false; end;
% %     if (final_pass), break; end; % This was the final pass and nothing has merged
% %     if (~something_merged), final_pass=true; end; % If we are done, do one last pass for final redistributes
% % end
% end %func


%--------------------------------------------------------------------------
% function centers = compute_centers_(X,labels)
% % JJJ optimized, isosplit
% [M,N]=size(X);
% centers=zeros(M,N);
% counts=accumarray(labels',1,[N,1])';
% for m=1:M
%     centers(m,:)=accumarray(labels',X(m,:)',[N,1])';
% end
% viLabel = find(counts);
% centers(:,viLabel)= bsxfun(@rdivide, centers(:,viLabel), counts(viLabel));
% end %func


%--------------------------------------------------------------------------
function export_diff_(P)
% export to _diff.bin, _diff.prb, _diff.prm files
error('not implemented yet');
if ~P.fTranspose_bin
    mnWav1 = reshape(load_bin_(P.vcFile, P.vcDataType), [], P.nChans);
    mnWav1 = mnWav1(:,P.viSite2Chan);
else
    mnWav1 = reshape(load_bin_(P.vcFile, P.vcDataType), P.nChans, []);
    mnWav1 = mnWav1(P.viSite2Chan, :)';
end

% mnWav1: nT x nSites

% fields to update, copy and save
P1 = P;
P1.vcFile = strrep(P.vcFile, '.bin', '_diff.bin');
P1.vcFile_prm = strrep(P.vcFile_prm, '.prm', '_diff.prm');
P1.probe_file = strrep(P.vcFile_prm, '.prm', '_diff.prb');
P1.fTranspose_bin = 0;
P1.vcCommonRef = 'none';
P1.fDetectBipolar = 1;
P1.nSites_ref = 0;

% differentiate channels and write to bin file (two column type)
nSites = numel(P.viSite2Chan);
viChan_HP = 1:2:nSites;
viChan_HN = 2:2:nSites;
viChan_VP = 1:(nSites-4);
viChan_VN = viChan_VP + 4;
viChan_P = [viChan_HP(1), toRow_([viChan_VP(1:2:end); viChan_HP(2:end-1); viChan_VP(2:2:end)]), viChan_HP(end)];
viChan_N = [viChan_HN(1), toRow_([viChan_VN(1:2:end); viChan_HN(2:end-1); viChan_VN(2:2:end)]), viChan_HN(end)];
mnWav2 = mnWav1(:,viChan_P) - mnWav1(:,viChan_N);
P1.nChans = size(mnWav2, 2);

% Output files
copyfile(P.vcFile_prm, P1.vcFile_prm, 'f');
edit_prm_file_(P1, P1.vcFile_prm);
write_bin_(P1.vcFile, mnWav2);
% write to probe file
% P1.probe_file
% mnWav2 = load_bin_(strrep(P.vcFile, '.bin', '_diff.bin'), P.vcDataType); %read back test


end %func


%--------------------------------------------------------------------------
function vr = toRow_(mr)
vr = mr(:)';
end


%--------------------------------------------------------------------------
function vr = toCol_(mr)
vr = mr(:);
end


%--------------------------------------------------------------------------
function vnWav1_mean = mean_excl_(mnWav1, P)
% calculate mean after excluding viSiteZero
viSiteZero = get_(P, 'viSiteZero');
if isempty(viSiteZero)
    vnWav1_mean = int16(mean(mnWav1,2));
else
    nSites_all = size(mnWav1, 2);
    nSites_excl = numel(viSiteZero);
    nSites = nSites_all - nSites_excl;
    vnWav1_mean = int16((sum(mnWav1,2) - sum(mnWav1(:,nSites_excl),2)) / nSites);
end
end %func


%--------------------------------------------------------------------------
function nCols = nColumns_probe_(P)
% Checkerboard four-column is considered as two column probe since
% two sites per vertical step
viShank_site = get_(P, 'viShank_site');
if ~isempty(viShank_site)
    viSites = find(P.viShank_site == P.viShank_site(1));    
    vrSiteY = P.mrSiteXY(viSites,2);
else
    vrSiteY = P.mrSiteXY(:,2);
end
vrSiteY_unique = unique(vrSiteY);
vnSites_group = hist(vrSiteY, vrSiteY_unique);
nCols = median(vnSites_group);
end %func


%--------------------------------------------------------------------------
function [mnWav2, cviSite_mean] = meanSite_drift_(mnWav1, P, viSite_repair)
% [mnWav1, cviSite_mean] = meanSite_drift_(mnWav1, P)
% [mnWav1, cviSite_mean] = meanSite_drift_(mnWav1, P, viSite_repair) %bad site repair

% mnWav1 = sites_repair_(mnWav1, P); % must repair site to perform vertical averaging
% this corrects for viSiteZero automatically
nSites = size(mnWav1,2);
nCols = nColumns_probe_(P);
viSiteZero = get_(P, 'viSiteZero');
cviSite_mean = cell(1, nSites);
viSites = 1:nSites;
fSingleShank = isSingleShank_(P);
if nargin<3
    viSite_repair = 1:nSites;
else
    viSite_repair = toRow_(viSite_repair);
end
mnWav2 = zeros(size(mnWav1), 'like', mnWav1);
for iSite = viSite_repair
    if fSingleShank
        viSite_shank1 = viSites; %faster        
    else
        viSite_shank1 = viSites(P.viShank_site == P.viShank_site(iSite));
    end
    for iNeigh=1:4
        viSite_mean1 = iSite + [-1,0,0,1] * nCols * iNeigh;
        viSite_mean1 = viSite_mean1(ismember(viSite_mean1, viSite_shank1));
        viSite_mean1(ismember(viSite_mean1, viSiteZero)) = [];
        if ~isempty(viSite_mean1), break; end
    end
    cviSite_mean{iSite} = viSite_mean1;
    mnWav2(:, iSite) = mean(mnWav1(:, viSite_mean1), 2);
end
end %func


%--------------------------------------------------------------------------
function flag = isSingleShank_(P)
viShank_site = get_(P, 'viShank_site');
if isempty(viShank_site)
    flag = 1;
else
    flag = numel(unique(viShank_site)) == 1;
end
end


%--------------------------------------------------------------------------
% function mrFet_clu = compute_centers_(S_clu)
% % calculate mean cluster waveform, using unique waveforms
% global tnWav_spk
% dimm_spk = size(tnWav_spk);
% mrFet_clu = zeros(dimm_spk(2), S_clu.nClu, 'single');
% % my feature is average cluster waveform
% for iClu=1:S_clu.nClu
%     mrFet_clu(:,iClu) = ;
% end %for
% end %func


%--------------------------------------------------------------------------
% function mrDist_clu = make_dists_matrix_(mrFet_clu, P)
% 
% end %func


%--------------------------------------------------------------------------
% function [inds1, inds2] = get_pairs_to_compare_(mrDist_clu)
% % find mutual min-dist pairs, remove inf
% 
% end %func


%--------------------------------------------------------------------------
% function [S_clu, centers, viPair_merged] = compare_pairs_(S_clu, centers, viClu1, viClu2)
% % find mutual min-dist pairs: correlation based merging only
% nPairs = numel(viClu1);
% for iPair = 1:nPairs
%     
% end 
% 
% centers = update_centers_(centers);
% end %func


%--------------------------------------------------------------------------
% function mrDist_clu = update_dists_matrix_(mrDist_clu, inds1, inds2, vi_inds_merged)
% % find mutual min-dist pairs
% mrDist_clu(sub2ind(size(mrDist_clu), inds1, inds2)) = inf; % compared dist ignored
% 
% end %func


%--------------------------------------------------------------------------
% function [inds1,inds2] = get_pairs_to_compare_(dists, thresh)
% [~,N]=size(dists); % square matrix, nClu
% inds1=[];
% inds2=[];
% % dists=make_dists_matrix(centers);
% % dists(find(comparisons_made(:)))=inf;
% % for j=1:N
% %     dists(j,j)=inf;
% % end;
% % important to only take the mutal closest pairs -- unlike how we originally did it
% %something_changed=1;
% %while (something_changed)
%     %something_changed=0;
%     [~,best_inds]=min(dists,[],1);
%     for j=1:N
%         if (best_inds(j)>j)
%             if (best_inds(best_inds(j))==j) % mutual
%                 if (dists(j,best_inds(j))<thresh)
%                     inds1(end+1)=j;
%                     inds2(end+1)=best_inds(j);
%                     dists(j,:)=inf;
%                     dists(:,j)=inf;
%                     dists(best_inds(j),:)=inf;
%                     dists(:,best_inds(j))=inf;
%                     %something_changed=1;
%                 end
%             end
%         end     
%     end
% %end;
% end %func


%--------------------------------------------------------------------------
function vhText = text_nClu_(S_clu, hAx)
nClu = numel(S_clu.vnSpk_clu);
viSite_clu = double(S_clu.viSite_clu);
y_offset = .3;
vhText = zeros(nClu+1,1);
for iClu=1:nClu
    n1 = S_clu.vnSpk_clu(iClu);
%         vcText1 = sprintf('%d, %0.1f', n1, n1/t_dur);
    vcText1 = sprintf('%d', n1); %show numbers
    vhText(iClu) = text(iClu, viSite_clu(iClu) + y_offset, vcText1, ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Parent', hAx);
end
vhText(end) = text(0,0,'#spk', 'VerticalAlignment', 'bottom', 'Parent', hAx);
end %func


%--------------------------------------------------------------------------
function mnWav2 = sgfilt_old_(mnWav, nDiff_filt, fGpu)
% works for a vector, matrix and tensor

if nargin<3, fGpu = isGpu_(mnWav); end
n1 = size(mnWav,1);
if n1==1, n1 = size(mnWav,2);  end
[via1, via2, via3, via4, vib1, vib2, vib3, vib4] = sgfilt4_(n1, fGpu);

if isvector(mnWav)
    switch nDiff_filt
        case 1
            mnWav2 = mnWav(via1) - mnWav(vib1);
        case 2
            mnWav2 = 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
        case 3
            mnWav2 = 3*(mnWav(via3) - mnWav(vib3)) + 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
        otherwise
            mnWav2 = 4*(mnWav(via4) - mnWav(vib4)) + 3*(mnWav(via3) - mnWav(vib3)) + 2*(mnWav(via2) - mnWav(vib2)) + mnWav(via1) - mnWav(vib1);
    end %switch
elseif ismatrix(mnWav)
    switch nDiff_filt
        case 1
            mnWav2 = mnWav(via1,:) - mnWav(vib1,:);
        case 2
            mnWav2 = 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
        case 3
            mnWav2 = 3*(mnWav(via3,:) - mnWav(vib3,:)) + 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
        otherwise
            mnWav2 = 4*(mnWav(via4,:) - mnWav(vib4,:)) + 3*(mnWav(via3,:) - mnWav(vib3,:)) + 2*(mnWav(via2,:) - mnWav(vib2,:)) + mnWav(via1,:) - mnWav(vib1,:);
    end %switch
else
    switch nDiff_filt
        case 1
            mnWav2 = mnWav(via1,:,:) - mnWav(vib1,:,:);
        case 2
            mnWav2 = 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
        case 3
            mnWav2 = 3*(mnWav(via3,:,:) - mnWav(vib3,:,:)) + 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
        otherwise
            mnWav2 = 4*(mnWav(via4,:,:) - mnWav(vib4,:,:)) + 3*(mnWav(via3,:,:) - mnWav(vib3,:,:)) + 2*(mnWav(via2,:,:) - mnWav(vib2,:,:)) + mnWav(via1,:,:) - mnWav(vib1,:,:);
    end %switch
end
end %func


%--------------------------------------------------------------------------
function mnWav1 = sgfilt_(mnWav, nDiff_filt, fGpu)
% works for a vector, matrix and tensor

if nargin<3, fGpu = isGpu_(mnWav); end
n1 = size(mnWav,1);
if n1==1, n1 = size(mnWav,2);  end
if nDiff_filt==0, mnWav1 = mnWav; return; end

[miA, miB] = sgfilt_init_(n1, nDiff_filt, fGpu);

if isvector(mnWav)
    mnWav1 = mnWav(miA(:,1)) - mnWav(miB(:,1));
    for i=2:nDiff_filt
        mnWav1 = mnWav1 + i * (mnWav(miA(:,i)) - mnWav(miB(:,i)));
    end
elseif ismatrix(mnWav)
    mnWav1 = mnWav(miA(:,1),:) - mnWav(miB(:,1),:);
    for i=2:nDiff_filt
        mnWav1 = mnWav1 + i * (mnWav(miA(:,i),:) - mnWav(miB(:,i),:));
    end
else
    mnWav1 = mnWav(miA(:,1),:,:) - mnWav(miB(:,1),:,:);
    for i=2:nDiff_filt
        mnWav1 = mnWav1 + i * (mnWav(miA(:,i),:,:) - mnWav(miB(:,i),:,:));
    end
end
end %func


%--------------------------------------------------------------------------
function flag = isGpu_(vr)
try
    flag = isa(vr, 'gpuArray');
catch
    flag = 0;
end
end


%--------------------------------------------------------------------------
function [tr, miRange] = mn2tn_gpu_(mr, spkLim, viTime, viSite)
% what to do if viTime goves out of the range?
% gpu memory efficient implementation
% it uses GPU if mr is in GPU

if nargin<4, viSite=[]; end %faster indexing
% if nargin<5, fMeanSubt=0; end

% JJJ 2015 Dec 24
% vr2mr2: quick version and doesn't kill index out of range
% assumes vi is within range and tolerates spkLim part of being outside
% works for any datatype
if isempty(viTime), tr=[]; return; end
[nT, nSites] = size(mr);
if ~isempty(viSite)
    mr = mr(:, viSite);
    nSites = numel(viSite); 
end
if iscolumn(viTime), viTime = viTime'; end

fGpu = isGpu_(mr);
viTime = gpuArray_(viTime, fGpu);
spkLim = gpuArray_(spkLim, fGpu);

viTime0 = [spkLim(1):spkLim(end)]'; %column
miRange = bsxfun(@plus, int32(viTime0), int32(viTime));
miRange = min(max(miRange, 1), nT);
% miRange = miRange(:);
tr = zeros([numel(viTime0), numel(viTime), nSites], 'int16');
dimm_tr = size(tr);
for iSite = 1:nSites
    if fGpu
%         vr1 = gpuArray(mr(:,iSite));
%         tr(:,:,iSite) = gather(vr1(miRange));
        tr(:,:,iSite) = gather_(reshape(mr(miRange, iSite), dimm_tr(1:2)));
    else
        tr(:,:,iSite) = reshape(mr(miRange, iSite), dimm_tr(1:2));
    end
end
end %func


%--------------------------------------------------------------------------
function [mr, fGpu] = gpuArray_(mr, fGpu)
if nargin<2, fGpu = 1; end
if ~fGpu, return; end
try
    mr = gpuArray(mr);
    fGpu = 1;
catch        
    try % retry after resetting the GPU memory
        gpuDevice(1); 
        mr = gpuArray(mr);
        fGpu = 1;
    catch % no GPU device found            
        fGpu = 0;
    end
end
end


%--------------------------------------------------------------------------
function nBytes_load = file_trim_(fid, nBytes_load, P)
if isempty(P.tlim_load) || ~P.fTranspose_bin, return; end

bytesPerSample = bytesPerSample_(P.vcDataType);
nSamples = floor(nBytes_load / bytesPerSample / P.nChans);

% Apply limit to the range of samples to load
nlim_load = min(max(round(P.tlim_load * P.sRateHz), 1), nSamples);
nSamples_load = diff(nlim_load) + 1;
nBytes_load = nSamples_load * bytesPerSample * P.nChans;
if nlim_load(1)>1, fseek_(fid, nlim_load(1), P); end
end %func


%--------------------------------------------------------------------------
function fseek_(fid_bin, iSample_bin, P)
% provide # of samples to skip for transpose multi-channel type
% Index starts with 1 unlike fseek() function

if ~P.fTranspose_bin, return; end
header_offset = get_set_(P, 'header_offset', 0);
iOffset = (iSample_bin-1) * P.nChans * bytesPerSample_(P.vcDataType) + header_offset;
if iOffset<0, iOffset = 0; end
fseek(fid_bin, iOffset, 'bof');
end %func


%--------------------------------------------------------------------------
function c = corr_vr_(vr1,vr2)
% Vectorize matrix and compute correlation
vr1 = vr1(:);
vr2 = vr2(:);
vr1 = vr1 - mean(vr1);
vr2 = vr2 - mean(vr2);
c = mean(vr1 .* vr2) / std(vr1,1) / std(vr2,1);
end %func


%--------------------------------------------------------------------------
% 10/8/17 JJJ: find correlation
function vrDist12 = xcorr2_mr_(mrWav1, mrWav2, arg1, arg2)
% vrDist12 = xcorr_mr_(mrWav1, mrWav2, nShift)
% vrDist12 = xcorr_mr_(mrWav1, mrWav2, cvi1, cvi2)
fMeanSubt_post = 1;
fSquared = 0;
if nargin == 3
    nShift = arg1;    
    nT = size(mrWav1, 1);
    [cvi1, cvi2] = shift_range_(nT, nShift);
else
    cvi1 = arg1;
    cvi2 = arg2;
end
if fSquared
    mrWav1 = mrWav1 .^ 2;
    mrWav2 = mrWav2 .^ 2;
end
vrDist12 = zeros(size(cvi1));
for iDist = 1:numel(vrDist12)    
    mr1 = mrWav1(cvi1{iDist},:);
    mr2 = mrWav2(cvi2{iDist},:);
    if ~fMeanSubt_post
        mr1 = meanSubt_(mr1);
        mr2 = meanSubt_(mr2);
    end
    vrDist12(iDist) = corr_(mr1(:), mr2(:), fMeanSubt_post);
end
end %func


%--------------------------------------------------------------------------
function [cvi1, cvi2] = shift_range_(nT, nShift, viShift)
% return ranges of two matrix to be time shifted
% [cvi1, cvi2] = shift_range_(nT, nShift): -nShift:nShift
% [cvi1, cvi2] = shift_range_(nT, vnShift)
if ~isempty(nShift)
    [cvi1, cvi2] = deal(cell(nShift*2+1, 1));
    viShift = -nShift:nShift;
else
    [cvi1, cvi2] = deal(cell(numel(viShift), 1));
end
viRange = 1:nT;
for iShift_ = 1:numel(viShift)
    iShift = viShift(iShift_);
    iShift1 = -round(iShift/2);
    iShift2 = iShift + iShift1;
    viRange1 = viRange + iShift1;
    viRange2 = viRange + iShift2;
    vl12 = (viRange1>=1 & viRange1<=nT) & (viRange2>=1 & viRange2<=nT);
    cvi1{iShift_} = viRange1(vl12);
    cvi2{iShift_} = viRange2(vl12);
end
end %func


%--------------------------------------------------------------------------
function S0 = save_log_(vcCmd, S0)

% save cluster info and save to file (append)
% check for crash
% todo: save differential increment from the beginning

if nargin<2, S0 = get0_(); end
[cS_log, P, S_clu, miClu_log] = get_(S0, 'cS_log', 'P', 'S_clu', 'miClu_log');

if ~isempty(strfind(vcCmd, 'annotate'))
    S_log = cS_log{end};
    S_log.csNote_clu = S_clu.csNote_clu;
    cS_log{end} = S_log;
else
    S_log = struct_('vcCmd', vcCmd, 'datenum', now(), 'csNote_clu', S_clu.csNote_clu);
    if isempty(cS_log) || ~iscell(cS_log)
        cS_log = {S_log};        
    else
        cS_log{end+1} = S_log;
    end
end

% Keep P.MAX_LOG history
if isempty(miClu_log)
    miClu_log = zeros([numel(S_clu.viClu), P.MAX_LOG], 'int16'); 
end
miClu_log(:, 2:end) = miClu_log(:, 1:end-1);
miClu_log(:, 1) = int16(S_clu.viClu);

%struct_save_(strrep(P.vcFile_prm, '.prm', '_log.mat'), 'cS_log', cS_log);
S_log.viClu = int16(S_clu.viClu);
struct_save_(S_log, strrep(P.vcFile_prm, '.prm', '_log.mat'), 0);
S0.cS_log = cS_log;
S0.miClu_log = miClu_log;
set(0, 'UserData', S0);

% update revert to list
ui_update_log_(cS_log);
end %func


%--------------------------------------------------------------------------
function S = struct_(varargin)
% smart about dealing with cell input
for iArg = 1:2:numel(varargin)
    try
        S.(varargin{iArg}) = varargin{iArg+1};
    catch
        disperr_('struct_');
    end
end
end


%--------------------------------------------------------------------------
function ui_update_log_(cS_log)
% the last one is selected
% persistent mh_history

% List recent activities
% if nargin<2, S0 = get0_(); end
% set(hMenu_history, 'Label', sprintf('Undo %s', cS_log.csCmd{end}), 'Enable', 'on');

% if isempty(mh_history)
mh_history = findobj('Type', 'uimenu', 'Tag', 'History');
P = get0_('P');
% Delete children and update
delete(mh_history.Children); %kill all children
for iMenu = 1:numel(cS_log) % reverse order
    iLog = numel(cS_log) - iMenu + 1;
    S_log1 = cS_log{iLog};
    vcLabel1 = sprintf('%s: %s', datestr(S_log1.datenum), S_log1.vcCmd);
    fEnable = (iMenu <= P.MAX_LOG) && iMenu~=1; 
    uimenu(mh_history, 'Label', vcLabel1, 'Callback', @(h,e)restore_log_(iMenu), ...
        'Checked', ifeq_(iMenu==1, 'on', 'off'), ...
        'Enable', ifeq_(fEnable, 'on', 'off'));
end
% update undo/redo menu
end %func


%--------------------------------------------------------------------------
function restore_log_(iMenu1)
% persistent mh_history
figure_wait_(1);
[cS_log, miClu_log, P] = get0_('cS_log', 'miClu_log', 'P');
S_clu1 = cS_log{end - iMenu1 + 1}; % last ones shown first
S_clu1.viClu = int32(miClu_log(:,iMenu1));

hMsg = msgbox_(sprintf('Restoring to %s (%s)', S_clu1.vcCmd, datestr(S_clu1.datenum)), 0);
[S_clu, S0] = S_clu_new_(S_clu1);

% Update checkbox
mh_history = findobj('Type', 'uimenu', 'Tag', 'History');
vhMenu = mh_history.Children;
vhMenu = vhMenu(end:-1:1); %reverse order
for iMenu = 1:numel(vhMenu)
    fTargetItem = iMenu==iMenu1;
    fEnable = ~fTargetItem && iMenu <= P.MAX_LOG;
    set(vhMenu(iMenu), ...
        'Checked', ifeq_(fTargetItem, 'on', 'off'), ...
        'Enable', ifeq_(fEnable, 'on', 'off'));
end %for

% update GUI
S0 = gui_update_(S0, S_clu);
% plot_FigWav_(S0); %redraw plot
% S0.iCluCopy = min(S0.iCluCopy, S_clu.nClu);
% S0.iCluPaste = [];
% set(0, 'UserData', S0);
% update_plot_(S0.hPaste, nan, nan); %remove paste cursor
% S0 = update_FigCor_(S0);        
% S0 = button_CluWav_simulate_(S0.iCluCopy, [], S0);
% set(0, 'UserData', S0);
close_(hMsg);
figure_wait_(0);
end %func


%--------------------------------------------------------------------------
function S0 = gui_update_(S0, S_clu)
if nargin<1, S0 = get0_(); end
if nargin<2, S_clu = S0.S_clu; end
plot_FigWav_(S0); %redraw plot
S0.iCluCopy = min(S0.iCluCopy, S_clu.nClu);

% S0.iCluCopy = 1;
S0.iCluPaste = [];
set(0, 'UserData', S0);
update_plot_(S0.hPaste, nan, nan); %remove paste cursor
S0 = update_FigCor_(S0);        
S0 = button_CluWav_simulate_(S0.iCluCopy, [], S0);
keyPressFcn_cell_(get_fig_cache_('FigWav'), 'z');
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
function vcFile_prm = makeprm_template_(vcFile_bin, vcFile_template, vcFile_prb)
% output prm file from a template file
% vcFile_prm is [vcFile_bin, vcFile_prb, '.prm'], after removing .bin and .prb extensions
csLines_prm = {};
csLines_prm{end+1} = sprintf('vcFile = ''%s'';', vcFile_bin);
csLines_prm{end+1} = sprintf('template_file = ''%s'';', vcFile_template);
if ~isempty(vcFile_prb)
    csLines_prm{end+1} = sprintf('probe_file = ''%s'';', vcFile_prb);
else
    S_prm = file2struct_(vcFile_template);
    vcFile_prb = S_prm.probe_file;
end

% update from meta file if exists
S_meta = read_meta_file_(subsFileExt_(vcFile_bin, '.meta'));
if ~isempty(S_meta)
    csLines_prm{end+1} = sprintf('nChans = %d;', S_meta.nChans);
    csLines_prm{end+1} = sprintf('sRateHz = %d;', S_meta.sRateHz);
    csLines_prm{end+1} = sprintf('uV_per_bit = %f;', S_meta.uV_per_bit);
end

% name prm file
[~,vcPostfix,~] = fileparts(vcFile_prb);
vcFile_prm = subsFileExt_(vcFile_bin, ['_', vcPostfix, '.prm']);

cellstr2file_(vcFile_prm, csLines_prm);
end


%--------------------------------------------------------------------------
function batch_mat_(vcFile_batch_mat, vcCommand)
% batch process binary file from a template file
% batch_(myfile_batch.mat, vcCommand)
%  file must contain: csFiles_bin, csFiles_template
%    optional: vrDatenum, datenum_start, csFiles_prb

if ~contains(lower(vcFile_batch_mat), '_batch.mat')
    fprintf(2, 'Must provide _batch.mat file format');
    return;
end
if nargin<2, vcCommand = ''; end
if isempty(vcCommand), vcCommand = 'spikesort'; end

% Input file format
S_batch = load(vcFile_batch_mat);
csFiles_bin = S_batch.csFiles_bin;
if ischar(csFiles_bin), csFiles_bin = {csFiles_bin}; end
nFiles = numel(csFiles_bin);
csFiles_template = get_(S_batch, 'csFiles_template');
if ischar(csFiles_template), csFiles_template = repmat({csFiles_template}, size(csFiles_bin)); end
csFiles_prb = get_(S_batch, 'csFiles_prb');
if isempty(csFiles_prb), csFiles_prb = ''; end
if ischar(csFiles_prb), csFiles_prb = repmat({csFiles_prb}, size(csFiles_bin)); end
csFiles_prm = cell(size(csFiles_bin));

for iFile = 1:nFiles
    try
        vcFile_prm1 = makeprm_template_(csFiles_bin{iFile}, csFiles_template{iFile}, csFiles_prb{iFile});
        fprintf('Created %s\n', vcFile_prm1);
        jrc3(vcCommand, vcFile_prm1);
        csFiles_prm{iFile} = vcFile_prm1;
    catch
        disperr_(sprintf('Failed to process %s', csFiles_bin{iFile}));
    end
end %for
S_batch.csFiles_prm = csFiles_prm;
write_struct_(vcFile_batch_mat, S_batch);
end %func


%--------------------------------------------------------------------------
function batch_plot_(vcFile_batch, vcCommand)
% vcFile_batch: .batch or _batch.mat file format (contains csFiles_prm)
% Collectively analyze multiple sessions
% error('not implemented yet');
if nargin<2, vcCommand=[]; end
if isempty(vcCommand), vcCommand='skip'; end %spikesort if doesn't exist

if ~exist(vcFile_batch, 'file'), fprintf(2, 'File does not exist\n'); return; end
if matchFileExt_(vcFile_batch, '.batch')
    edit(vcFile_batch); %show script
    csFiles_prm = importdata(vcFile_batch);
    % Removing comments that starts with "%"
    func_comment = @(vc)vc(1) == '%';
    viComment = cellfun(@(vc)func_comment(strtrim(vc)), csFiles_prm);
    csFiles_prm(viComment) = [];
end


% run the sorting and collect data. quantify the quality
cS_plot_file = cell(size(csFiles_prm));
for iFile=1:numel(csFiles_prm)
    try
        vcFile_prm1 = csFiles_prm{iFile};
        jrc3('clear');
        if ~strcmpi(vcCommand, 'skip')                        
            jrc3(vcCommand, vcFile_prm1);
            S0 = get0_();
        else
            S0 = load_cached_(vcFile_prm1);
        end
        cS_plot_file{iFile} = S_plot_new_(S0);
    catch
        disp(lasterr());
    end
end %for

% plot cS_plot_file
S_plot_show(cS_plot_file); % save to _batch.mat (?)

end %func


%--------------------------------------------------------------------------
function S_plot = S_plot_new_(S0)
% S_plot contains quantities to be plotted
% Copied from jrclust.m quality_metric_
global tnWav_spk
if nargin<1, S0 = get0_(); end
P = S0.P;

vrVrms_site = single(S0.vrThresh_site(:)) / P.qqFactor;
vrSnr_evt = single(abs(S0.vrAmp_spk(:))) ./ vrVrms_site(S0.viSite_spk(:));
t_dur = double(max(S0.viTime_spk) - min(S0.viTime_spk)) / P.sRateHz;
vrRate_site = cellfun(@numel, S0.cviSpk_site)' / t_dur;
% nSites = numel(S0.cviSpk_site);

% calc # spikes exceeding detection threshold
vnSite_evt = zeros(size(S0.viTime_spk), 'int16');
for iSite = 1:numel(S0.cviSpk_site)
    viSpk_site1 = S0.cviSpk_site{iSite};
    mrMin_site1 = squeeze_(min(tnWav_spk(:,:,viSpk_site1)));
    vrThresh_site1 = -abs(S0.vrThresh_site(P.miSites(:, iSite)));
    vnSite_evt(viSpk_site1) = sum(bsxfun(@lt, mrMin_site1, vrThresh_site1(:)));
end

% cluster meta analysis (cluster of clusters)


% Compute cluster stats
mrMin_clu = uV2bit_(squeeze_(min(S0.S_clu.trWav_spk_clu)));
vrSnr_clu = abs(mrMin_clu(1,:))' ./ vrVrms_site(S0.S_clu.viSite_clu);
vrRate_clu = cellfun(@numel, S0.S_clu.cviSpk_clu)' / t_dur;
mrThresh_clu = -abs(S0.vrThresh_site(P.miSites(:,S0.S_clu.viSite_clu)));
vnSite_clu = sum(mrMin_clu < mrThresh_clu)';

S_plot = makeStruct_(vrVrms_site, vrRate_site, t_dur, P, ...
    vrSnr_evt, vnSite_evt, vrSnr_clu, vrRate_clu, vnSite_clu);
end %func


%--------------------------------------------------------------------------
function save_var_(vcFile, varargin)
% must pass 
struct_save_(struct_(varargin{:}), vcFile);
end


%--------------------------------------------------------------------------
function S_plot_show(cS_plot_file)
% plot clusters of clusters
% @TODO
end %func


%--------------------------------------------------------------------------
function [S_clu, nClu_merged] = S_clu_pv_merge_(S_clu, P) %update mrWavCor when you merge
global tnWav_spk tnWav_raw
MAD_THRESH = -4;
nClu_merged = 0;
[viSite_spk] = get0_('viSite_spk');
% mrWavCor = S_clu.mrWavCor;
nClu = S_clu.nClu;
fprintf('S_clu_pv_merge_\n');

% Identify clusters to remove, update and same (no change), disjoint sets
[vrMinDist_logz_clu, viMinDist_clu] = S_clu_pca_dist_(S_clu);
%vi_clu1 = find(vrMinDist_logz_clu < MAD_THRESH);
vi_clu1 = find(vrMinDist_logz_clu < -max(vrMinDist_logz_clu));
if isempty(vi_clu1), return; end

vi_clu2 = viMinDist_clu(vi_clu1);
viMap_clu = 1:nClu;
viMap_clu(vi_clu1) = vi_clu2;
viClu_same = setdiff(1:nClu, union(vi_clu1, vi_clu2));
viClu_remove = setdiff(1:nClu, viMap_clu);
viClu_update = setdiff(setdiff(1:nClu, viClu_same), viClu_remove);
% viClu_update = setdiff(1:nClu, viClu_same);

% update cluster number
try S_clu.icl(viClu_remove) = []; catch, end
S_clu = S_clu_map_index_(S_clu, viMap_clu); %index mapped
P.fVerbose = 0;
S_clu = S_clu_refrac_(S_clu, P); % remove refrac spikes

% update cluster waveforms and distance
S_clu = clu2wav_(S_clu, viSite_spk, tnWav_spk, tnWav_raw, viClu_update); %update cluster waveforms
S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update);
% S_clu = S_clu_refresh_(S_clu); % remove empty and remap
S_clu = S_clu_remove_empty_(S_clu);

nClu_merged = nClu - S_clu.nClu;
fprintf('\n\tnClu: %d->%d (%d merged)\n', nClu, S_clu.nClu, nClu_merged);
end %func


%--------------------------------------------------------------------------
function [vrMinDist_logz_clu, viMinDist_clu] = S_clu_pca_dist_(S_clu)
global tnWav_raw tnWav_spk

MAX_REAL_DIST = 50;
MAX_SAMPLE = 2000;
fUseMean = 1; %use median instead
nPc = 2;
fUseRaw = 1;
fUsePvCorr = 1;
nShift = 6;
fUseSd = 1;

P = S_clu.P;

trWav_clu = ifeq_(fUseRaw, S_clu.trWav_raw_clu, S_clu.trWav_spk_clu);
if ~fUseMean
    viSite_spk = get0_('viSite_spk');
end
nClu = S_clu.nClu;
nSamples = size(trWav_clu,1);
% mrPv1_clu = zeros(nSamples, nClu);
% mrPv1_clu = zeros(nSamples, nClu);
nDelay = 3;
[mrPv1_clu, mrPv2_clu, mrPv3_clu] = deal(zeros(size(trWav_clu,1), nClu));
for iClu=1:nClu
%     [~, mrPv1_clu(:,iClu)] = pca(trWav_clu(:,:,iClu), 'NumComponents', 1);
    if fUseMean
        mrWav_clu1 = trWav_clu(:,:,iClu);        
    else
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        viSpk_clu1 = viSpk_clu1(viSite_spk(viSpk_clu1) == S_clu.viSite_clu(iClu));
        viSpk_clu1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
        if fUseRaw
%             mrWav_clu1 = single(median(tnWav_raw(:,:,viSpk_clu1), 3));
            mrWav_clu1 = single(reshape(tnWav_raw(:,:,viSpk_clu1), nSamples, []));
        else
%             mrWav_clu1 = single(median(tnWav_spk(:,:,viSpk_clu1), 3));
            mrWav_clu1 = single(reshape(tnWav_spk(:,:,viSpk_clu1), nSamples, []));
        end
    end
    if fUseSd                              
        mrPv1_clu(:,iClu) = std(mrWav_clu1,1,2);
        if nPc>=2, mrPv2_clu(:,iClu) = mr_std2_(mrWav_clu1, nDelay)'; end
        if nPc>=2, mrPv3_clu(:,iClu) = mr_std2_(mrWav_clu1, nDelay*2)'; end
    else
        [~, mrPv_clu1] = pca(mrWav_clu1, 'NumComponents', nPc);
        mrPv1_clu(:,iClu) = mrPv_clu1(:,1);
        if nPc>=2, mrPv2_clu(:,iClu) = mrPv_clu1(:,2); end
        if nPc>=3, mrPv3_clu(:,iClu) = mrPv_clu1(:,3); end
    end
end

if fUsePvCorr
    func1 = @(x)1 - max(abs(xcorr_mr_(x, nShift)),[],3);
    switch nPc
        case 1
            mrPcDist_clu = func1(mrPv1_clu);
        case 2
            mrPcDist_clu = (func1(mrPv1_clu) + func1(mrPv2_clu))/2;
        case 3
            mrPcDist_clu = (func1(mrPv1_clu) + func1(mrPv2_clu) + func1(mrPv3_clu))/3;
    end          
else
    [vrPc1_clu, vrPv1_clu] = pca(mrPv1_clu, 'NumComponents', 1);
    switch nPc
        case 1            
            mrPc_clu = vrPc1_clu;
        case 2
            [vrPc2_clu, vrPv2_clu] = pca(mrPv2_clu, 'NumComponents', 1);
            mrPc_clu = [vrPc1_clu, vrPc2_clu];
        case 3
            [vrPc2_clu, vrPv2_clu] = pca(mrPv2_clu, 'NumComponents', 1);
            [vrPc3_clu, vrPv3_clu] = pca(mrPv3_clu, 'NumComponents', 1);
            mrPc_clu = [vrPc1_clu, vrPc2_clu, vrPc3_clu];
    end    
    mrPcDist_clu = pdist2_(abs(mrPc_clu));
end
mrRealDist_clu = pdist2_(P.mrSiteXY(S_clu.viSite_clu,:));

%mrPcDist_clu(sub2ind([nClu,nClu], 1:nClu, 1:nClu)) = nan;
mrPcDist_clu(tril(true(nClu)) | mrRealDist_clu > MAX_REAL_DIST) = nan; %ignore bottom half

% lower triangle only
[vrMinDist_clu, viMinDist_clu] = min(mrPcDist_clu);
vrRealDist_clu = mrRealDist_clu(sub2ind([nClu,nClu], 1:nClu, viMinDist_clu));

vrMinDist_logz_clu = zeros(size(vrMinDist_clu));
vi_NotNan = find(vrMinDist_clu > 0);
vrMinDist_logz_clu(vi_NotNan) = zscore_(log(vrMinDist_clu(vi_NotNan)));
% vrMinDist_logz_clu(2:end) = madscore_(log(vrMinDist_clu(2:end)));

if nargout==0
    figure; plot(vrMinDist_logz_clu, vrRealDist_clu, '.'); 
    xlabel('min clu dist (log-MAD pc)'); ylabel('real clu dist (um)'); grid on;    
    vi_clu1 = find(vrMinDist_logz_clu < MAD_THRESH);
    vi_clu2 = viMinDist_clu(vi_clu1);
    vr_dist12 = vrMinDist_logz_clu(vi_clu1);
    arrayfun(@(a,b,c)fprintf('(%d,%d,%0.2f), ', a,b,c), vi_clu1, vi_clu2, vr_dist12);
    fprintf('\n');
end
end %func


%--------------------------------------------------------------------------
function [sd2, viRange1, viRange2] = mr_std2_(mr, nDelay)
if nDelay==0, sd2 = std(mr,1,2); return; end

% determine shift
nT = size(mr,1);
iShift1 = -round(nDelay/2);
iShift2 = nDelay + iShift1;
viRange1 = max((1:nT) + iShift1, 1);
viRange2 = min((1:nT) + iShift2, nT);

% viRange1 = (1:nT) + iShift1;
% viRange2 = (1:nT) + iShift2;
% vl12 = (viRange1>=1 & viRange1<=nT) & (viRange2>=1 & viRange2<=nT);
% viRange1 = viRange1(vl12);
% viRange2 = viRange2(vl12);

mr1 = mr(viRange1,:);
mr2 = mr(viRange2,:);

sd2 = sqrt(abs(mean(mr1.*mr2,2) - mean(mr1,2).*mean(mr2,2)));  
end %func


%--------------------------------------------------------------------------
function sd2 = mr_std3_(mr, viRange1, viRange2)
mr1 = mr(viRange1,:);
mr2 = mr(viRange2,:);
sd2 = sqrt(abs(mean(mr1.*mr2,2) - mean(mr1,2).*mean(mr2,2)));  
end %func


%--------------------------------------------------------------------------
function trCorr = xcorr_mr_(mrPv_clu, nShift)
% vrDist12 = xcorr_mr_(mrWav1, mrWav2, nShift)
% vrDist12 = xcorr_mr_(mrWav1, mrWav2, cvi1, cvi2)
if nShift==0
    trCorr = corr_(mrPv_clu);
    return;
end
nT = size(mrPv_clu,1);
[cvi1, cvi2] = shift_range_(nT, nShift);
nClu = size(mrPv_clu,2);
trCorr = zeros([nClu, nClu, numel(cvi1)], 'like', mrPv_clu);
for iShift = 1:numel(cvi1)    
    trCorr(:,:,iShift) = corr_(mrPv_clu(cvi1{iShift},:), mrPv_clu(cvi2{iShift},:));
end
end %func


%--------------------------------------------------------------------------
function mr = madscore_(mr)
% maximum absolute difference transformation

mr = bsxfun(@minus, mr, median(mr));
vr = median(abs(mr));
mr = bsxfun(@rdivide, mr, vr);
end %func


%--------------------------------------------------------------------------
function [mrPc1, mrPc2] = pca_tr_(tn)
% returns first principal component across sites
persistent tn_
global tnWav_spk

if nargin<1, tn = tn_; end
if isempty(tn)
    tn_ = tnWav_spk;
    tn = tn_;
end

% mrPv = zeros(size(tr,1), size(tr,3), 'single');
mrPc1 = zeros(size(tn,2), size(tn,3), 'single');
mrPc2 = ifeq_(nargout>1, mrPc1, []);
        
% tr = single(tr);
% n = size(tn,2);
tr = meanSubt_tr_(single(tn));
nSpk = size(tn,3);
% tic
if isempty(mrPc2)
    parfor iSpk = 1:nSpk
    %     [~, mrPv(:,iSpk)] = pca(tr(:,:,iSpk), 'NumComponents', 1);
    %     mr1 = single(tn(:,:,iSpk));
    %     mr1 = bsxfun(@minus, mr1, mean(mr1));
        mr1 = tr(:,:,iSpk);
        [V,D] = eig(mr1*mr1');
        D = sqrt(diag(D));
        mrPc1(:,iSpk) = V(:,end)' * mr1 / D(end); %V(:,end)' * mr1;
    %     mrPc(:,iSpk) = mr1' * V(:,end) / sqrt(D(end)); %V(:,end)' * mr1;
    %     mrPc(:,iSpk) = pca(single(tn(:,:,iSpk)), 'NumComponents', 1); %equivalent
    end
else
    parfor iSpk = 1:nSpk
        mr1 = tr(:,:,iSpk);
        [V,D] = eig(mr1*mr1');
        D = sqrt(diag(D));
        mrPc1(:,iSpk) = V(:,end)' * mr1 / D(end); %V(:,end)' * mr1;
        mrPc2(:,iSpk) = V(:,end-1)' * mr1 / D(end-1); %V(:,end)' * mr1;
    end
end
% toc
end %func


%--------------------------------------------------------------------------
function mrPc1 = pc1_tr_(tn)
% returns first principal component across sites

mr0 = single(squeeze_(tn(:,1,:)));
vrPv0 = zscore_(pca(mr0', 'NumComponents', 1));
dimm_tn = size(tn);
mr = single(reshape(tn, dimm_tn(1), []));
mr = bsxfun(@minus, mr, mean(mr));
mrPc1 = reshape(vrPv0' * mr, dimm_tn(2:3));
end %func


%--------------------------------------------------------------------------
function tr = meanSubt_tr_(tr)
dimm = size(tr);
mr = reshape(tr,size(tr,1),[]);
mr = bsxfun(@minus, mr, mean(mr));
tr = reshape(mr, dimm);
end %func


%--------------------------------------------------------------------------
function export_(varargin)
% export_(): export S0 struct to the workspace
% export_(var1, var2): export fields in S0 struct to the workspace
nArgs = nargin();
csVars = varargin;
csVars = csVars(~isempty_(csVars));
if isempty(csVars), csVars = {'S0'}; end
S0 = get0_();
for iArg = 1:numel(csVars)
    vcVar = csVars{iArg};
    if isempty(vcVar), continue; end
    if ~strcmpi(vcVar, 'S0')
        var = get_(S0, vcVar);
    else
        var = S0;
    end
    if isempty(var)
        fprintf(2, '''%s'' does not exist\n', vcVar);
    else
        assignin('base', vcVar, var);
        fprintf('assigned ''%s'' to workspace\n', vcVar);
    end
end
end %func


%--------------------------------------------------------------------------
function vl = isempty_(cvr)
if iscell(cvr)
    vl = cellfun(@isempty, cvr);
else
    vl = isempty(cvr);
end
end %func


%--------------------------------------------------------------------------
function S0 = clear_log_(S0)
S0.cS_log = {};
S0.miClu_log = [];
set0_(S0);
delete_files_(strrep(S0.P.vcFile_prm, '.prm', '_log.mat'), 0);
end %func


%--------------------------------------------------------------------------
function auto_(P)
% S0 = get0_();

S0 = load_cached_(P); % load cached data or from file if exists
S_clu = get_(S0, 'S_clu');
S_clu.P = P;

if isempty(S_clu)
    fprintf(2, 'You must sort first by running "jrc3 sort".\n');
    return;
end
[S_clu, S0] = post_merge_(S_clu, P);
S0 = clear_log_(S0);
save0_(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
end %func


%--------------------------------------------------------------------------
function ydB = pow2db_(y)
ydB = (10.*log10(y)+300)-300;
end %func


%--------------------------------------------------------------------------
function mr12 = pdist2_(mr1, mr2)
% mr12 = pdist2_(mr1) % self distance
% mr12 = pdist2_(mr1, mr2)
% mr1: n1xd, mr2: n2xd, mr12: n1xn2

% mr12 = sqrt(eucl2_dist_(mr1', mr2'));
% 20% faster than pdist2 for 10000x10 x 100000x10 single
if nargin==2
    mr12 = sqrt(bsxfun(@plus, sum(mr2'.^2), bsxfun(@minus, sum(mr1'.^2)', 2*mr1*mr2')));
else
    vr1 = sum(mr1'.^2);
    mr12 = sqrt(bsxfun(@plus, vr1, bsxfun(@minus, vr1', 2*mr1*mr1')));
end
end %func


%--------------------------------------------------------------------------
function z = zscore_(x, flag, dim)
if isempty(x), z=[]; return; end
if nargin < 2, flag = 0; end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);
end %func


%--------------------------------------------------------------------------
function C = corr_(A, B, fMeanSubt)
% mr = corr_(A, B)
% mr = corr_(A) % n1 x n2 becomes n1 x 
% mr = corr_(vr1, vr2) % single coefficient
if nargin<3, fMeanSubt = 1; end
% https://stackoverflow.com/questions/9262933/what-is-a-fast-way-to-compute-column-by-column-correlation-in-matlab
if fMeanSubt, A = bsxfun(@minus,A,mean(A)); end %% zero-mean
A = bsxfun(@times,A,1./sqrt(sum(A.^2))); %% L2-normalization
if nargin == 1
    C = A' * A;
else
    if fMeanSubt, B = bsxfun(@minus,B,mean(B)); end %%% zero-mean
    B = bsxfun(@times,B,1./sqrt(sum(B.^2))); %% L2-normalization
    C = A' * B;
end
end %func


%--------------------------------------------------------------------------
% 6/23/JJJ

%--------------------------------------------------------------------------
function traces_(P, fDebug_ui_, vcFileId)
% show raw traces
% If file format is nChans x nSamples, load subset of file (fTranspose=1)
% If file format is nSamples x nChans, load all and save to global (fTranspose=0)
% 2017/6/22 James Jun: Added multiview (nTime_traces )
global fDebug_ui
global mnWav mnWav1 % only use if P.fTranspose=0
if nargin==0
    P = get0_('P'); 
else
    set0_(P);
end
if nargin<2, fDebug_ui_=0; end
if nargin<3, vcFileId=''; end
if isempty(P), disperr_('traces_: P is empty'); return; end
% S0 = load0_(P);
fDebug_ui = fDebug_ui_;
S0 = load_cached_(P, 0);
set(0, 'UserData', S0);
set0_(fDebug_ui);
% S0 = load_cached_(P, 0); %don't load raw waveform

% get file to show
iFile_show = 1; %files to display for clustered together
if ~isempty(P.csFile_merge)
    csFiles_bin = filter_files_(P.csFile_merge);
    if numel(csFiles_bin)==1
        vcFile_bin = csFiles_bin{1}; 
    else %show multiple files        
        if isempty(vcFileId)
            arrayfun(@(i)fprintf('%d: %s\n', i, csFiles_bin{i}), 1:numel(csFiles_bin), 'UniformOutput', 0);
            fprintf('---------------------------------------------\n');
            vcFileId = input('Please specify Fild ID from the list above:', 's');
        end
        if isempty(vcFileId), return; end
        iFile_show = str2num(vcFileId);        
        try
            vcFile_bin = csFiles_bin{iFile_show};
        catch
            return;
        end
    end
else
    vcFile_bin = P.vcFile; % if multiple files exist, load first
end
set0_(iFile_show);
tlim_bin = P.tlim;
% if isempty(nTime_traces), nTime_traces = 1; end

% Open file
fprintf('Opening %s\n', vcFile_bin);
[fid_bin, nBytes_bin] = fopen_(vcFile_bin, 'r');
if isempty(fid_bin), fprintf(2, '.bin file does not exist: %s\n', vcFile_bin); return; end
nSamples_bin = floor(nBytes_bin / bytesPerSample_(P.vcDataType) / P.nChans);
nLoad_bin = min(round(diff(tlim_bin) * P.sRateHz), nSamples_bin);
if tlim_bin(1)>0
    iSample_bin = ceil(tlim_bin(1) * P.sRateHz) + 1; %offset sample number    
else
    iSample_bin = 1; %sample start location
end    
nlim_bin = [0,nLoad_bin-1] + iSample_bin;
if nlim_bin(1) < 1, nlim_bin = [1, nLoad_bin]; end
if nlim_bin(2) > nSamples_bin, nlim_bin = [-nLoad_bin+1, 0] + nSamples_bin; end

nTime_traces = get_(P, 'nTime_traces');
[cvn_lim_bin, viRange_bin] = sample_skip_(nlim_bin, nSamples_bin, nTime_traces);    
if P.fTranspose_bin   
    mnWav = [];
    fseek_(fid_bin, iSample_bin, P);
    if nTime_traces > 1            
        mnWav1 = load_bin_multi_(fid_bin, cvn_lim_bin, P)';        
    else        
        mnWav1 = load_bin_(fid_bin, P.vcDataType, [P.nChans, nLoad_bin])'; %next keypress: update tlim_show    
    end
%     @TODO: load from cvn_lim_bin specifiers. check for end or beginning when keyboard command
else %load whole thing
    mnWav = load_bin_(fid_bin, P.vcDataType, [nSamples_bin, P.nChans]); %next keypress: update tlim_show
    fclose(fid_bin);
    fid_bin = []; 
    %mnWav1 = mnWav((nlim_bin(1):nlim_bin(2)), :);
    mnWav1 = mnWav(viRange_bin, :);
    disp('Entire raw traces are cached to RAM since fTranspose=0.');
end %if

% full screen width
hFig_traces = create_figure_('Fig_traces', [0 0 .5 1], vcFile_bin); %remove all other figure traces
hAx = axes_new_(hFig_traces); % create axis
hPlot = line(hAx, nan, nan, 'Color', [1 1 1]*.5, 'Parent', hAx, 'LineWidth', .5);
hPlot_edges = plot(nan, nan, 'Color', [1 0 0]*.5, 'Parent', hAx, 'LineWidth', 1);
set(hAx, 'Position',[.05 .05 .9 .9], 'XLimMode', 'manual', 'YLimMode', 'manual');
S_fig = makeStruct_(hAx, hPlot, nlim_bin, fid_bin, nSamples_bin, nLoad_bin, hPlot_edges);
S_fig.maxAmp = P.maxAmp;
S_fig.vcTitle = '[H]elp; (Sft)[Up/Down]:Scale(%0.1f uV); (Sft)[Left/Right]:Time; [F]ilter; [J]ump T; [C]han. query; [R]eset view; [P]SD; [S]pike; [A]ux chan; [E]xport; [T]race; [G]rid';
S_fig.csHelp = { ...    
    'Left/Right: change time (Shift: x4)', ...
    '[J]ump T', ...
    '[Home/End]: go to beginning/end of file', ...
    '---------', ...
    'Up/Down: change scale (Shift: x4)', ...
    'Zoom: Mouse wheel', ...
    '[x/y/ESC]: zoom direction', ...
    'Pan: hold down the wheel and drag', ...
    '[R]eset view', ...
    '---------', ...
    '[F]ilter toggle', ...
    '[S]pike toggle', ...
    'Gri[D] toggle', ...
    '[T]races toggle', ... 
    '---------', ...    
    '[C]hannel query', ...
    '[A]ux channel display', ...
    '[P]ower spectrum', ...    
    '[E]xport to workspace', ...
    };
S_fig = struct_append_(S_fig, ...
    struct('vcGrid', 'on', 'vcFilter', 'off', 'vcSpikes', 'on', 'vcTraces', 'on'));
set(hFig_traces, 'UserData', S_fig);
set(hFig_traces, 'color', 'w', 'KeyPressFcn', @keyPressFcn_Fig_traces_, 'BusyAction', 'cancel', 'CloseRequestFcn', @close_hFig_traces_);
mouse_figure(hFig_traces);

Fig_traces_plot_(1); % Plot spikes and color clusters 
end %func


%--------------------------------------------------------------------------
% plot data
function Fig_traces_plot_(fAxis_reset)
% fAxis_reset: reset the axis limit
% [usage]
% Fig_traces_plot_()
% Fig_traces_plot_(fAxis_reset)
% 2017/06/22 James Jun
% 6/22 JJJ: added seperator lines, fixed the reset view and spike view

global mnWav1 mrWav1 % current timeslice to plot
if nargin<1, fAxis_reset = 0; end
fWait = msgbox_('Plotting...',0,1);
fShuttleOrder = 1; %shuffle cluster color
[S0, P, S_clu] = get0_();
[hFig, S_fig] = get_fig_cache_('Fig_traces'); 
figure_wait_(1, hFig); drawnow;
sRateHz = P.sRateHz / P.nSkip_show;
viSamples1 = 1:P.nSkip_show:size(mnWav1,1);
spkLim = round(P.spkLim / P.nSkip_show); %show 2x of range
P.vcFilter = get_filter_(P);
if strcmpi(S_fig.vcFilter, 'on')
    P1=P; P1.sRateHz = sRateHz; P1.fGpu = 0;
    P1.vcFilter = get_set_(P, 'vcFilter_show', P.vcFilter);
    if P.fft_thresh>0, mnWav1 = int16(fft_clean_(mnWav1, P)); end
    mrWav1 = bit2uV_(filt_car_(mnWav1(viSamples1, P.viSite2Chan), P1), P1);
    vcFilter_show = P1.vcFilter;
else
    mrWav1 = meanSubt_(single(mnWav1(viSamples1, P.viSite2Chan))) * P.uV_per_bit;    
    vcFilter_show = 'off';
end
viSites = 1:numel(P.viSite2Chan);
% mrWav1 = meanSubt_(single(mnWav1(:, P.viSite2Chan))) * P.uV_per_bit;
% hide bad channels
nTime_traces = get_(P, 'nTime_traces');
if isempty(nTime_traces) || nTime_traces==1
    vrTime_bin = ((S_fig.nlim_bin(1):P.nSkip_show:S_fig.nlim_bin(end))-1) / P.sRateHz;    
    vcXLabel = 'Time (s)';
else    
    vrTime_bin = (0:(size(mrWav1,1)-1)) / (P.sRateHz / P.nSkip_show) + (S_fig.nlim_bin(1)-1) / P.sRateHz;
    [cvn_lim_bin, viRange_bin, viEdges] = sample_skip_(S_fig.nlim_bin, S_fig.nSamples_bin, nTime_traces);
    tlim_show = (cellfun(@(x)x(1), cvn_lim_bin([1,end]))) / P.sRateHz;
    vcXLabel = sprintf('Time (s), %d segments merged (%0.1f ~ %0.1f s, %0.2f s each)', nTime_traces, tlim_show, diff(P.tlim));
    mrX_edges = vrTime_bin(repmat(viEdges(:)', [3,1]));
    mrY_edges = repmat([0;numel(P.viSite2Chan)+1;nan],1,numel(viEdges));
    set(S_fig.hPlot_edges, 'XData', mrX_edges(:), 'YData', mrY_edges(:));
    csTime_bin = cellfun(@(x)sprintf('%0.1f', x(1)/P.sRateHz), cvn_lim_bin, 'UniformOutput', 0);
    set(S_fig.hAx, {'XTick', 'XTickLabel'}, {vrTime_bin(viEdges), csTime_bin});
end
multiplot(S_fig.hPlot, S_fig.maxAmp, vrTime_bin, mrWav1, viSites);
% axis(S_fig.hAx, [vrTime_bin(1), vrTime_bin(end), viSites(1)-1, viSites(end)+1]);
grid(S_fig.hAx, S_fig.vcGrid);
set(S_fig.hAx, 'YTick', viSites);
title_(S_fig.hAx, sprintf(S_fig.vcTitle, S_fig.maxAmp));
xlabel(S_fig.hAx, vcXLabel); 
ylabel(S_fig.hAx, 'Site #'); 
set(S_fig.hPlot, 'Visible', S_fig.vcTraces);

% Delete spikes from other threads
S_fig_ = get(hFig, 'UserData');
if isfield(S_fig_, 'chSpk'), delete_multi_(S_fig_.chSpk); end
if isfield(S_fig, 'chSpk'), delete_multi_(S_fig.chSpk); end

% plot spikes
if strcmpi(S_fig.vcSpikes, 'on') && isfield(S0, 'viTime_spk')
    viTime_spk = S0.viTime_spk - int32(S0.viT_offset_file(S0.iFile_show));
    if nTime_traces > 1
        viSpk1 = find(in_range_(viTime_spk, cvn_lim_bin));
        [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
        viTime_spk1 = round(reverse_lookup_(viTime_spk1, viRange_bin) / P.nSkip_show);
    else
        viSpk1 = find(viTime_spk >= S_fig.nlim_bin(1) & viTime_spk < S_fig.nlim_bin(end));
        [viSite_spk1, viTime_spk1] = multifun_(@(vr)vr(viSpk1), S0.viSite_spk, viTime_spk);
        viTime_spk1 = round((viTime_spk1 - S_fig.nlim_bin(1) + 1) / P.nSkip_show); %time offset
    end        
    t_start1 = single(S_fig.nlim_bin(1) - 1) / P.sRateHz;
    viSite_spk1 = single(viSite_spk1);
    % check if clustered
    if isempty(S_clu)
        nSites = size(mrWav1,2);
        chSpk = cell(nSites, 1);
        for iSite=1:nSites %deal with subsample factor
            viSpk11 = find(viSite_spk1 == iSite);
            if isempty(viSpk11), continue; end
            viTime_spk11 = viTime_spk1(viSpk11);
            [mrY11, mrX11] = vr2mr3_(mrWav1(:,iSite), viTime_spk11, spkLim); %display purpose x2
%             vr2mr_spk_(mrWav1(:,iSite), viTime_spk11, P);
            mrT11 = single(mrX11-1) / sRateHz + t_start1;
            chSpk{iSite} = line(nan, nan, 'Color', [1 0 0], 'LineWidth', 1.5, 'Parent', S_fig.hAx);
            multiplot(chSpk{iSite}, S_fig.maxAmp, mrT11, mrY11, iSite);
        end        
    else % different color for each clu
        viClu_spk1 = S_clu.viClu(viSpk1);        
        mrColor_clu = [jet(S_clu.nClu); 0 0 0];        
        vrLineWidth_clu = (mod((1:S_clu.nClu)-1, 3)+1)'/2 + .5;  %(randi(3, S_clu.nClu, 1)+1)/2;
        if fShuttleOrder
            mrColor_clu = shuffle_static_(mrColor_clu, 1);
            vrLineWidth_clu = shuffle_static_(vrLineWidth_clu, 1);
        end
        nSpk1 = numel(viTime_spk1);
        chSpk = cell(nSpk1, 1);
        for iSpk1 = 1:nSpk1
            iTime_spk11 = viTime_spk1(iSpk1);
            iSite11 = viSite_spk1(iSpk1);
            [mrY11, mrX11] = vr2mr3_(mrWav1(:,iSite11), iTime_spk11, spkLim); %display purpose x2
            mrT11 = double(mrX11-1) / sRateHz + t_start1;
            iClu11 = viClu_spk1(iSpk1);
            if iClu11<=0, continue; end
%                 vrColor1 = [0 0 0]; 
%                 lineWidth1 = .5;
%             else
            vrColor1 = mrColor_clu(iClu11,:);
            lineWidth1 = vrLineWidth_clu(iClu11);
%             end
            chSpk{iSpk1} = line(nan, nan, 'Color', vrColor1, 'LineWidth', lineWidth1, 'Parent', S_fig.hAx);
            multiplot(chSpk{iSpk1}, S_fig.maxAmp, mrT11, mrY11, iSite11);
        end
    end
    S_fig.chSpk = chSpk;
else
    % delete spikes    
    S_fig.chSpk = [];    
end
if fAxis_reset, fig_traces_reset_(S_fig); end
set(hFig, 'UserData', S_fig, 'Name', sprintf('%s: filter: %s', P.vcFile_prm, (vcFilter_show)));
figure_wait_(0, hFig);
close_(fWait);
end %func


%--------------------------------------------------------------------------
function fig_traces_reset_(S_fig)
global mnWav1

if nargin<1, [hFig, S_fig] = get_fig_cache_('Fig_traces');  end
% axis(S_fig.hAx, [S_fig.nlim_bin / P.sRateHz, 0, nSites+1]);
P = get0_('P');
nTime_traces = get_(P, 'nTime_traces');
if nTime_traces > 1
    tlim1 = ([0, size(mnWav1,1)] + S_fig.nlim_bin(1) - 1) / P.sRateHz;
    tlim1 = round(tlim1*1000)/1000;
    axis(S_fig.hAx, [tlim1, 0, numel(P.viSite2Chan)+1]);
else
    axis(S_fig.hAx, [S_fig.nlim_bin / P.sRateHz, 0, numel(P.viSite2Chan)+1]);
end 
end %func


%--------------------------------------------------------------------------
function [cvnlim_bin, viRange, viEdges] = sample_skip_(nlim_bin, nSamples_bin, nTime_traces)
% return a limit that 
% nlim_bin=[81 90]; nSamples_bin=100; nTime_traces=5;
% edges to set to nan
% 2017/6/22 James Jun: Added nTime_traces multiview
% 6/23 JJJ: edge samples to be set to nan (gap)

if nTime_traces==1 || isempty(nTime_traces)
    cvnlim_bin = {nlim_bin};
    viRange = nlim_bin(1):nlim_bin(end);
    viEdges = [];
    return;
end
nSkip = floor(nSamples_bin / nTime_traces);
cvnlim_bin = arrayfun(@(i)nlim_bin + (i-1)*nSkip, 1:nTime_traces, 'UniformOutput', 0);
% modulus algebra wrap around
for i=1:nTime_traces
    lim1 = mod(cvnlim_bin{i}-1, nSamples_bin)+1;
    if lim1(1) > lim1(2)
        lim1 = [1, diff(nlim_bin)+1];
    end
    cvnlim_bin{i} = lim1;
end
if nargout>=2
    viRange = cell2mat_(cellfun(@(x)x(1):x(2), cvnlim_bin, 'UniformOutput', 0));
end
if nargout>=3 %compute the number of samples
    viEdges = cumsum(cellfun(@(x)diff(x)+1, cvnlim_bin));
    viEdges = [1, viEdges(1:end-1)];
%     viEdges = sort([viEdges, viEdges+1], 'ascend'); %two sample gaps
end
end %func


%--------------------------------------------------------------------------
function vl = in_range_(vi, cvi)
vl = false(size(vi));
if ~iscell(cvi), cvi = {cvi}; end
for i=1:numel(cvi)
    lim1 = cvi{i};
    vl = vl | (vi >= lim1(1) & vi <= lim1(2));
end
end %func


%--------------------------------------------------------------------------
function [viA, vl] = reverse_lookup_(viB, viA2B)
% viB must belong to viA2B
% viB=[3 1 1 5 3 3 2]; viA2B=[1 3 5];

viB = int32(viB);
viA2B = int32(viA2B);
[vl, viA] = ismember(int32(viB), int32(viA2B));

% vl = ismember(viB, viA2B);
% assert(all(vl), 'reverse_lookup_: all viB must belong to viA2B');
% viA = arrayfun(@(i)find(viA2B==i), viB(vl), 'UniformOutput', 1); 
end %func


%--------------------------------------------------------------------------
function [tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk, viTime_spk, vnAmp_spk, vnThresh_site] = ...
    wav2spk_(mnWav1, vrWav_mean1, P, viTime_spk, viSite_spk, mnWav1_pre, mnWav1_post)
% tnWav_spk: spike waveform. nSamples x nSites x nSpikes
% trFet_spk: nSites x nSpk x nFet
% miSite_spk: nSpk x nFet
% spikes are ordered in time
% viSite_spk and viTime_spk is uint32 format, and tnWav_spk: single format
% mnWav1: raw waveform (unfiltered)
% wav2spk_(mnWav1, vrWav_mean1, P)
% wav2spk_(mnWav1, vrWav_mean1, P, viTime_spk, viSite_spk)
% 6/27/17 JJJ: accurate spike detection at the overlap region
% 6/29/17 JJJ: matched filter supported

if nargin<4, viTime_spk = []; end
if nargin<5, viSite_spk = []; end
if nargin<6, mnWav1_pre = []; end
if nargin<7, mnWav1_post = []; end
[tnWav_spk_raw, tnWav_spk, trFet_spk, miSite_spk] = deal([]);

nFet_use = get_set_(P, 'nFet_use', 2);

fMerge_spk = 1; %debug purpose
fShift_pos = 0; % shift center position based on center of mass
fRecenter_spk = 0;

% Filter
fprintf('\tFiltering spikes...'); t_filter = tic;
if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
    mnWav1 = [mnWav1_pre; mnWav1; mnWav1_post];
end
if P.fft_thresh>0
    [mnWav1, P.fGpu] = fft_clean_(mnWav1, P); 
    mnWav1 = int16(mnWav1);
end
[mnWav2, vnWav11] = filt_car_(mnWav1, P);

% common mode rejection
if P.blank_thresh > 0
    if isempty(vnWav11)
        vnWav11 = mr2ref_(mnWav2, P.vcCommonRef, P.viSiteZero); %vrWav_mean1(:);    
    end
    vlKeep_ref = car_reject_(vnWav11(:), P);
    fprintf('Rejecting %0.3f %% of time due to motion\n', (1-mean(vlKeep_ref))*100 );
else
    vlKeep_ref = [];
end
% set0_(vlKeep_ref);
fprintf('\ttook %0.1fs\n', toc(t_filter));

if isempty(get_(P, 'vcFilter_detect'))
    mnWav3 = mnWav2;
else
    [mnWav3, nShift_post] = filter_detect_(mnWav2, P); % apply stdfilter for peak search
end

% detect spikes or use the one passed from the input (importing)
nPad_pre = size(mnWav1_pre,1);
vnThresh_site = get0_('vnThresh_site');
if isempty(vnThresh_site)
    vnThresh_site = gather_(int16(mr2rms_(mnWav3, 1e5) * P.qqFactor));
end
if isempty(viTime_spk) || isempty(viSite_spk)
    P_ = setfield(P, 'nPad_pre', nPad_pre);
    [viTime_spk, vnAmp_spk, viSite_spk] = detect_spikes_(mnWav3, vnThresh_site, vlKeep_ref, P_);
else
    viTime_spk = viTime_spk + nPad_pre;
    vnAmp_spk = mnWav3(sub2ind(size(mnWav3), viTime_spk, viSite_spk)); % @TODO read spikes at the site and time
end
vnAmp_spk = gather_(vnAmp_spk);
% if nShift_post~=0, viTime_spk = viTime_spk + nShift_post; end % apply possible shift due to filtering

% reject spikes within the overlap region
if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)    
    ilim_spk = [nPad_pre+1, size(mnWav3,1) - size(mnWav1_post,1)]; %inclusive
    viKeep_spk = find(viTime_spk >= ilim_spk(1) & viTime_spk <= ilim_spk(2));
    [viTime_spk, vnAmp_spk, viSite_spk] = multifun_(@(x)x(viKeep_spk), viTime_spk, vnAmp_spk, viSite_spk);    
end%if
if isempty(viTime_spk), return; end


% Extract spike waveforms and build a spike table
fprintf('\tExtracting features'); t_fet = tic;
mnWav2 = gather_(mnWav2); %do in CPU. 10.2s in GPU, 10.4s in CPU
if fRecenter_spk % center site is where the energy is the highest, if disabled min is chosen
    tnWav_spk = mn2tn_wav_spk2_(mnWav2, viSite_spk, viTime_spk, P);
    %[~, viMaxSite_spk] = max(squeeze_(std(single(tnWav_spk))));
    [~, viMaxSite_spk] = max(squeeze_(max(tnWav_spk) - min(tnWav_spk)));
    viSite_spk = P.miSites(sub2ind(size(P.miSites), viMaxSite_spk(:), viSite_spk));
end
[tnWav_spk_raw, tnWav_spk] = mn2tn_wav_(mnWav1, mnWav2, viSite_spk, viTime_spk, P); fprintf('.');
nSite_use = P.maxSite*2+1 - P.nSites_ref;
if nSite_use==1, nFet_use=1; end
assert(nSite_use >0, 'nSites_use = maxSite*2+1 - nSites_ref must be greater than 0');
switch nFet_use
    case 3
        [viSite2_spk, viSite3_spk] = find_site_spk23_(tnWav_spk, viSite_spk, P); fprintf('.');
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        mrFet2 = trWav2fet_(mn2tn_wav_spk2_(mnWav2, viSite2_spk, viTime_spk, P), P); fprintf('.');
        mrFet3 = trWav2fet_(mn2tn_wav_spk2_(mnWav2, viSite3_spk, viTime_spk, P), P); fprintf('.');
        trFet_spk = permute(cat(3, mrFet1, mrFet2, mrFet3), [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [viSite_spk(:), viSite2_spk(:), viSite3_spk(:)]; %nSpk x nFet
    case 2
        viSite2_spk = find_site_spk23_(tnWav_spk, viSite_spk, P); fprintf('.');
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        mrFet2 = trWav2fet_(mn2tn_wav_spk2_(mnWav2, viSite2_spk, viTime_spk, P), P); fprintf('.');
        trFet_spk = permute(cat(3, mrFet1, mrFet2), [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [viSite_spk(:), viSite2_spk(:)]; %nSpk x nFet
    case 1
        mrFet1 = trWav2fet_(tnWav_spk, P); fprintf('.');
        trFet_spk = permute(mrFet1, [1,3,2]); %nSite x nFet x nSpk
        miSite_spk = [viSite_spk(:)];
    otherwise
        error('wav2spk_: nFet_use must be 1, 2 or 3');
end
if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end

fprintf('\ttook %0.1fs\n', toc(t_fet));
end %func


%--------------------------------------------------------------------------
function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_fast_(vrWav1, P, thresh1)
% P: spkThresh, qqSample, qqFactor, fGpu, uV_per_bit
% vrWav1 can be either single or int16
% 6/27/17 JJJ: bugfix: hard set threshold is applied

% Determine threshold
MAX_SAMPLE_QQ = 300000; 
fSpikeRefrac_site = 0;
if nargin < 2, thresh1 = []; end
if ~isempty(P.spkThresh), thresh1 = int16(P.spkThresh); end
if thresh1==0, [viSpk1, vrSpk1] = deal([]); return; end % bad site
if isempty(thresh1)
%     if P.nDiff_filt>0 %already centered
        thresh1 = median(abs(subsample_vr_(vrWav1, MAX_SAMPLE_QQ)));
        thresh1 = int16(single(thresh1)* P.qqFactor / 0.6745);
%     else %uncentered
%         vrWav1_sub = subsample_vr_(vrWav1, MAX_SAMPLE_QQ);
%         med1 = median(vrWav1_sub);
%         thresh1 = median(abs(vrWav1_sub - med1));
%         thresh1 = int16(single(thresh1)* P.qqFactor / 0.6745) + abs(med1);
%     end
end
% end
thresh1 = int16(thresh1);

% detect valley turning point. cannot detect bipolar
% pick spikes crossing at least three samples
nneigh_min = get_(P, 'nneigh_min_detect');
viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min);
if P.fDetectBipolar
   viSpk1 = [viSpk1; find_peak_(-vrWav1, thresh1, nneigh_min)]; 
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

% apply spike merging on the same site
nRefrac = int32(abs(P.spkRefrac));
if P.refrac_factor > 1
    nRefrac = int32(round(double(nRefrac) * P.refrac_factor));
end
if isGpu_(viSpk1)
    [viSpk1, vrSpk1, thresh1] = multifun_(@gather, viSpk1, vrSpk1, thresh1);
end
if fSpikeRefrac_site %perform spike refractive period per site (affects exact mode)
    [viSpk1, vrSpk1] = spike_refrac_(viSpk1, vrSpk1, [], nRefrac); %same site spikes
end
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_position_(S_clu, viClu_update)
% determine cluster position from spike position
% 6/27/17 JJJ: multiple features supported (single dimension such as energy and Vpp)

global trFet_spk
if nargin<2, viClu_update = []; end
P = S_clu.P;
if ~isfield(S_clu, 'vrPosX_clu'), S_clu.vrPosX_clu = []; end
if ~isfield(S_clu, 'vrPosY_clu'), S_clu.vrPosY_clu = []; end

if isempty(S_clu.vrPosX_clu) || ~isempty(S_clu.vrPosY_clu)
    viClu_update = [];
end
if isempty(viClu_update)
    [vrPosX_clu, vrPosY_clu] = deal(zeros(S_clu.nClu, 1));
    viClu1 = 1:S_clu.nClu;
else % selective update
    vrPosX_clu = S_clu.vrPosX_clu;
    vrPosY_clu = S_clu.vrPosY_clu;
    viClu1 = viClu_update(:)';
end
viSites_fet = 1:(1+P.maxSite*2-P.nSites_ref);
for iClu = viClu1
%     viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    [viSpk_clu1, viSites_clu1] = S_clu_subsample_spk_(S_clu, iClu);
    if isempty(viSpk_clu1), continue; end
    
    viSites_clu1 = viSites_clu1(1:end-P.nSites_ref);
    mrVp1 = squeeze_(trFet_spk(viSites_fet,1,viSpk_clu1));
    mrSiteXY1 = single(P.mrSiteXY(viSites_clu1,:)); %electrode
    
    vrPosX_clu(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,1), 2));
    vrPosY_clu(iClu) = median(centroid_mr_(mrVp1, mrSiteXY1(:,2), 2));
end
S_clu.vrPosX_clu = vrPosX_clu;
S_clu.vrPosY_clu = vrPosY_clu;
end %func


%--------------------------------------------------------------------------
function [mn1, nShift_post] = filter_detect_(mn, P, vcMode)
% returns spatial sd
% mn0 = single(mn);
% mn0 = bsxfun(@minus, mn0, mean(mn0, 2)) .^ 2;
% 6/29/17 JJJ: filter detection

if nargin<3, vcMode = get_(P, 'vcFilter_detect'); end
% global vrFilt_spk nShift_post_

viSites_use = 1:(1+2*P.maxSite - P.nSites_ref);
viSites_ref = (1+2*P.maxSite - P.nSites_ref+1):(1+2*P.maxSite);
fprintf('filter_detect\n\t'); t1= tic;
miSites = gpuArray_(P.miSites(viSites_use, :));
miSites_ref = gpuArray_(P.miSites(viSites_ref, :));
nShift_post = 0;
switch lower(vcMode)
    case 'chancor'
        mn1 = chancor_(mn, P);        
    case 'matched'
        vrFilt_spk = get0_('vrFilt_spk');
        if isempty(vrFilt_spk)
            lim_ = round([3,5]/8 * size(mn,1));
            mn_ = mn(lim_(1):lim_(2),:);
            [vrFilt_spk, vrVaf, nShift_post] = calc_matched_filt_(mn_, P); %detect primary 
            set0_(vrFilt_spk, nShift_post);
        else
            nShift_post = get0_('nShift_post');
        end
%         nShift_post = nShift_post_;
        mn1 = int16(conv2(single(gather_(mn)), vrFilt_spk(:), 'same'));
%         mn1 = shift_mr_(mn1, nShift_post); % or do it later in the spike detection phase
%         figure; plot(xcorr(mn(:,41), mn1(:,41), 10));
%         nShift_post = P.spkLim(1)-1;
    case 'autocov'
        mn1 = -int16(filt_corr(single(mn), 2));
    case 'std-chan'
        mn1 = zeros(size(mn), 'like', mn);
        for iSite=1:size(mn,2)
            %mn_ = mn(:, P.miSites(viSites_use, iSite));
            %vn_ref = mean(mn_(:,viSites_ref),2);
            %mn_ = bsxfun(@minus, mn_(:,viSites_use), vn_ref);
           % mn1(:, iSite) = -int16(std(single(mn(:, P.miSites(:, iSite))), 0, 2));
            mn1(:, iSite) = -int16(std(single(mn(:, miSites(:,iSite))), 1, 2));
            %mn1(:,iSite) = mean(mn_.^2,2) - mean(mn_,2).^2;
            fprintf('.');
        end
        
    case 'std-time'
        % envelop filter. this affects the threshold. 
        
    case 'nmean'
        mn1 = zeros(size(mn), 'like', mn);
        for iSite=1:size(mn,2)
            %mn_ = mn(:, P.miSites(viSites_use, iSite));
            %vn_ref = mean(mn_(:,viSites_ref),2);
            %mn_ = bsxfun(@minus, mn_(:,viSites_use), vn_ref);
           % mn1(:, iSite) = -int16(std(single(mn(:, P.miSites(:, iSite))), 0, 2));
            %mn1(:, iSite) = mn(:,iSite) - int16(mean(single(mn(:,miSites_ref(:,iSite))),2));
            mn1(:, iSite) = int16(mean(mn(:,miSites(:,iSite)), 2));
            %mn1(:,iSite) = mean(mn_.^2,2) - mean(mn_,2).^2;
            fprintf('.');
        end
        
end %switch
fprintf('\n\ttook %0.1fs\n', toc(t1));
% mn1 = -int16(sqrt(mn1 / size(P.miSites,1)));
end %func


%--------------------------------------------------------------------------
function [tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, viTime_spk0, viSite_spk0)
% file loading routine. keep spike waveform (tnWav_spk) in memory
% assume that the file is chan x time format
% usage:
% [tnWav_raw, tnWav_spk, S0] = file2spk_(P)
%   
% [tnWav_raw, tnWav_spk, S0] = file2spk_(P, viTime_spk, viSite_spk)
%   construct spike waveforms from previous time markers
% 6/29/17 JJJ: Added support for the matched filter

% global vrFilt_spk; % mrPv_global
if nargin<2, viTime_spk0 = []; end
if nargin<3, viSite_spk0 = []; end
[tnWav_raw, tnWav_spk, trFet_spk, S0] = deal([]);

viTime_spk0 = viTime_spk0(:);
viSite_spk0 = viSite_spk0(:);

if isempty(P.csFile_merge)
    csFile = {P.vcFile};
else
    csFile = filter_files_(P.csFile_merge);
end
if ~isempty(get_(P, 'vcFile_thresh'))
    try
        S_thresh = load(P.vcFile_thresh);
        vnThresh_site = S_thresh.vnThresh_site;
        set0_(vnThresh_site);
        fprintf('Loaded %s\n', P.vcFile_thresh);
    catch
        disperr_('vcFile_thresh load error');
    end
end
if isempty(csFile), error('No binary files found.'); end
[tnWav_raw, tnWav_spk, trFet_spk, miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site] = deal({});    
viT_offset_file = zeros(size(csFile));
nFiles = numel(csFile);    
[nSamples1, nLoads] = deal(0); % initialize the counter
[vrFilt_spk, mrPv_global] = deal([]); % reset the template
set0_(mrPv_global, vrFilt_spk); % reeset mrPv_global and force it to recompute
for iFile=1:nFiles
    fprintf('File %d/%d: detecting spikes from %s\n', iFile, nFiles, csFile{iFile});
    t1 = tic;
    [fid1, nBytes_file1] = fopen_(csFile{iFile}, 'r');
    nBytes_file1 = file_trim_(fid1, nBytes_file1, P);
    [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file1, P);
%         nSamples1 = 0; %accumulated sample offset        
    viT_offset_file(iFile) = nSamples1;
    mnWav11_pre = [];
    for iLoad1 = 1:nLoad1
        fprintf('Processing %d/%d of file %d/%d...\n', iLoad1, nLoad1, iFile, nFiles);
        nSamples11 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);
        fprintf('\tLoading from file...'); t_load_ = tic;
        [mnWav11, vrWav_mean11] = load_file_(fid1, nSamples11, P);
        fprintf('took %0.1fs\n', toc(t_load_));
        if iLoad1 < nLoad1
            mnWav11_post = load_file_preview_(fid1, P);
        else
            mnWav11_post = [];
        end
        [viTime_spk11, viSite_spk11] = filter_spikes_(viTime_spk0, viSite_spk0, nSamples1 + [1, nSamples11]);
        [tnWav_raw{end+1}, tnWav_spk{end+1}, trFet_spk{end+1}, miSite_spk{end+1}, ...
            viTime_spk{end+1}, vrAmp_spk{end+1}, vnThresh_site{end+1}] ...
                = wav2spk_(mnWav11, vrWav_mean11, P, viTime_spk11, viSite_spk11, mnWav11_pre, mnWav11_post);
        viTime_spk{end} = viTime_spk{end} + nSamples1;
        nSamples1 = nSamples1 + nSamples11;
        if iLoad1 < nLoad1, mnWav11_pre = mnWav11(end-P.nPad_filt+1:end, :); end
        clear mnWav11 vrWav_mean11;
        nLoads = nLoads + 1;
    end %for
    fclose(fid1);   
    t_dur1 = toc(t1);
    t_rec1 = (nBytes_file1 / bytesPerSample_(P.vcDataType) / P.nChans) / P.sRateHz;
    fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
        iFile, nFiles, ...
        t_dur1, nBytes_file1/1e6, nBytes_file1/t_dur1/1e6, t_rec1/t_dur1);
end %for
% catch
%     disperr_();
% end

% combine files and 
tnWav_raw = cat(3, tnWav_raw{:});
tnWav_spk = cat(3, tnWav_spk{:});
trFet_spk = cat(3, trFet_spk{:});
[miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site] = ...
    multifun_(@(x)cat(1, x{:}), miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site);
vrThresh_site = mean(single(vnThresh_site),1);
viSite_spk = miSite_spk(:,1);
if size(miSite_spk,2) >= 2
    viSite2_spk = miSite_spk(:,2);
else
    viSite2_spk = [];
end

% set S0
dimm_raw = size(tnWav_raw); %2x longer spkLim_raw = spkLim*2;
dimm_spk = size(tnWav_spk);
dimm_fet = size(trFet_spk);
nSites = numel(P.viSite2Chan);
cviSpk_site = arrayfun(@(iSite)find(miSite_spk(:,1) == iSite), 1:nSites, 'UniformOutput', 0);
if size(miSite_spk,2) >= 2
    cviSpk2_site = arrayfun(@(iSite)find(miSite_spk(:,2) == iSite), 1:nSites, 'UniformOutput', 0);
else
    cviSpk2_site = cell(1, nSites);
end
if size(miSite_spk,2) >= 3
    cviSpk3_site = arrayfun(@(iSite)find(miSite_spk(:,3) == iSite), 1:nSites, 'UniformOutput', 0);
else
    cviSpk3_site = [];
end
S0 = makeStruct_(P, viSite_spk, viSite2_spk, viTime_spk, vrAmp_spk, vrThresh_site, dimm_spk, ...
    cviSpk_site, cviSpk2_site, cviSpk3_site, dimm_raw, viT_offset_file, dimm_fet, nLoads, ...
    mrPv_global, vrFilt_spk);
end %func


%--------------------------------------------------------------------------
function mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update)
% symmetric matrix and common basis comparison only
% 6/29/17 JJJ: distance-based neighboring cluster selection
fUsePeak2 = 1;
fWaveform_raw = 1;
fMaxSite_excl = 0; %excl max site that undergoes most change during drift
fDiff_raw = 0; %worse off if enabled higher FP/FN
nShift = 6; % +/-n number of samples to compare time shift
nSite_overlap_thresh = floor(P.maxSite);
maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
% maxDist_site_um = maxDist_site_um/2; % half the distance

if nargin<3, viClu_update = []; end
if ~isfield(S_clu, 'mrWavCor'), viClu_update = []; end
fprintf('Computing waveform correlation...\n\t'); t1 = tic;
tmrWav_clu = ifeq_(fWaveform_raw, S_clu.tmrWav_raw_clu, S_clu.tmrWav_spk_clu);
if fDiff_raw && fWaveform_raw, tmrWav_clu = sgfilt_(tmrWav_clu, P.nDiff_filt); end
maxSite = ceil(P.maxSite);
nClu = S_clu.nClu;
if fUsePeak2
    [viSite_clu, viSite2_clu, viSite3_clu] = S_clu_peak2_(S_clu);
else
    viSite_clu = S_clu.viSite_clu;
end
mrWavCor = zeros(nClu);
miSites = P.miSites;
% miSites = P.miSites(1:(maxSite+1), :); % center

nT = size(tmrWav_clu, 1);
[cviShift1, cviShift2] = shift_range_(nT, nShift);
viLags = 1:numel(cviShift1);
if isempty(viClu_update)
    vlClu_update = true(nClu, 1); 
else
    vlClu_update = false(nClu, 1);
    vlClu_update(viClu_update) = 1;
    mrWavCor0 = S_clu.mrWavCor;
    nClu_pre = size(mrWavCor0, 1);
    vlClu_update((1:nClu) > nClu_pre) = 1;
end
try
for iClu2 = 1:nClu         
    iSite_clu2 = viSite_clu(iClu2);    
    if iSite_clu2==0 || isnan(iSite_clu2), continue; end
    viSite2 = miSites(:,iSite_clu2); 
    if fMaxSite_excl, viSite2 = viSite2(2:end); end
    
    % find neighboring clu within physical radius. todo: second min site filtering
    if ~fUsePeak2
        viClu1 = find(ismember(viSite_clu, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));
    else
        viClu1 = find(viSite_clu == iSite_clu2 | viSite2_clu == iSite_clu2 | viSite3_clu == iSite_clu2 | ...
                viSite_clu == viSite2_clu(iClu2) | viSite_clu == viSite3_clu(iClu2)); %viSite2_clu == viSite2_clu(iClu2) |
    end
    viClu1(viClu1 <= iClu2) = []; % symmetric matrix comparison
    tmrWav_clu21 = tmrWav_clu(:,viSite2,:); %temp
    mrWav_clu21 = tmrWav_clu21(:,:,iClu2);
    tmrWav_clu21 = tmrWav_clu21(:,:,viClu1);
    for iClu11 = 1:numel(viClu1)
        iClu1 = viClu1(iClu11);
        if ~vlClu_update(iClu1) && ~vlClu_update(iClu2)
            mrWavCor(iClu1, iClu2) = mrWavCor0(iClu1, iClu2);
        else            
            iSite_clu1 = viSite_clu(iClu1);
            if iSite_clu1==0 || isnan(iSite_clu1), continue; end
            if iSite_clu1 == iSite_clu2
                mrWav_clu1 = tmrWav_clu21(:,:,iClu11);
                mrWav_clu2 = mrWav_clu21;
            else
                viSite1 = miSites(:,iSite_clu1);
                viSite12 = find(ismember(viSite2, viSite1));
                if isempty(viSite12), continue; end
%                 if numel(viSite12) < nSite_overlap_thresh, continue; end
                mrWav_clu1 = tmrWav_clu21(:,viSite12,iClu11);
                mrWav_clu2 = mrWav_clu21(:,viSite12);
            end
            mrWavCor(iClu1, iClu2) = max(xcorr2_mr_(mrWav_clu1, mrWav_clu2, cviShift1, cviShift2));
        end
    end
    fprintf('.');
end
catch
    disperr_();
end
mrWavCor = mrWavCor + mrWavCor'; %make it symmetric
mrWavCor(mrWavCor==0) = nan;
fprintf('\n\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [vrFilt_spk, vrVaf, nShift_post] = calc_matched_filt_(mnWav1, P) %detect primary 
% generate a matched filter Kernel
% determine the treshold
% 6/29/17 JJJ: Spike waveform matched fitler determination
% 6/30/17 JJJ: Range optimization

vnThresh_site = gather_(int16(mr2rms_(mnWav1, 1e5) * P.qqFactor));
[viTime_spk, vnAmp_spk, viSite_spk] = detect_spikes_(mnWav1, vnThresh_site, [], P);

% extract wave forms
nSpks = numel(viSite_spk);
nSites = numel(P.viSite2Chan);
%spkLim = [-1, 1] * round(mean(abs(P.spkLim))); %balanced
% spkLim = [-1, 1] * round(max(abs(P.spkLim))); %balanced
spkLim = [-1, 1] * round(min(abs(P.spkLim))); %balanced
spkLim = spkLim + [1,0]; 

mnWav_spk = zeros(diff(spkLim) + 1, nSpks, 'int16');
for iSite = 1:nSites
    viiSpk11 = find(viSite_spk == iSite);
    if isempty(viiSpk11), continue; end
    viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
    mnWav_spk(:,viiSpk11) = gather_(vr2mr3_(mnWav1(:,iSite), viTime_spk11, spkLim));
end

[vrFilt_spk, ~, vrVaf] = pca(single(mnWav_spk'), 'NumComponents', 1, 'Centered', 0);
vrFilt_spk = flipud(vrFilt_spk(:));
if abs(min(vrFilt_spk)) > abs(max(vrFilt_spk)), vrFilt_spk = -vrFilt_spk; end
% vrFilt_spk(1) = []; % start from 1 correction
vrFilt_spk = vrFilt_spk - mean(vrFilt_spk);

% vrFilt_spk = vrFilt_spk / (vrFilt_spk.'*vrFilt_spk);
vrVaf = cumsum(vrVaf);
vrVaf = vrVaf / vrVaf(end);
[~,nShift_post] = max(vrFilt_spk);
nShift_post = round(numel(vrFilt_spk)/2 - nShift_post);
end %func


%--------------------------------------------------------------------------
function [viSite, viSite2, viSite3] = S_clu_peak2_(S_clu)
mrMin_clu = squeeze_(min(S_clu.tmrWav_spk_clu) - max(S_clu.tmrWav_spk_clu));
% mrMin_clu = squeeze_(min(S_clu.tmrWav_spk_clu));
% mrMin_clu = squeeze_(min(S_clu.tmrWav_raw_clu));

[~, viSite] = min(mrMin_clu);
if nargout>=2
    mrMin_clu(sub2ind(size(mrMin_clu), viSite, 1:numel(viSite))) = 0;
    [~, viSite2] = min(mrMin_clu);
end
if nargout>=3
    mrMin_clu(sub2ind(size(mrMin_clu), viSite2, 1:numel(viSite2))) = 0;
    [~, viSite3] = min(mrMin_clu);
end
end %func


%--------------------------------------------------------------------------
function [mrFet1_clu1, iSite_clu1] = S_clu_getFet_(S_clu, iClu, viSite2_spk)
global trFet_spk
if nargin<3, viSite2_spk = get0_('viSite2_spk'); end
iSite_clu1 = S_clu.viSite_clu(iClu);
viSpk_clu1 = S_clu.cviSpk_clu{iClu};
if isempty(viSpk_clu1), mrFet1_clu1=[]; return; end
mrFet1_clu1 = squeeze_(trFet_spk(:,1,viSpk_clu1));
if size(trFet_spk,2) >= 2
    mrFet2_clu1 = squeeze_(trFet_spk(:,2,viSpk_clu1));
    viSpk_clu1_site2 = find(viSite2_spk(viSpk_clu1) == iSite_clu1);
    mrFet1_clu1(:,viSpk_clu1_site2) = mrFet2_clu1(:,viSpk_clu1_site2);
end
end %func
        

%--------------------------------------------------------------------------
function S_clu = S_clu_cleanup_(S_clu, P)
% 17/7/3: Cluster cleanup routine, Mahal distance based outlier removal
% global trFet_spk

viSite2_spk = get0_('viSite2_spk');
thresh_mad_clu = get_set_(P, 'thresh_mad_clu', 7.5);
if thresh_mad_clu == 0, return; end % aborted

fprintf('Cleaning up clusters\n\t'); t1=tic;
% viClu = S_clu.viClu;
warning off;
for iClu = 1:S_clu.nClu
    try
        mrFet1_clu1 = S_clu_getFet_(S_clu, iClu, viSite2_spk);
%         iSite_clu1 = S_clu.viSite_clu(iClu);
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
%         if isempty(viSpk_clu1), continue; end
%         mrFet1_clu1 = squeeze_(trFet_spk(:,1,viSpk_clu1));
%         mrFet2_clu1 = squeeze_(trFet_spk(:,2,viSpk_clu1));
%         viSpk_clu1_site2 = find(viSite2_spk(viSpk_clu1) == iSite_clu1);
%         mrFet1_clu1(:,viSpk_clu1_site2) = mrFet2_clu1(:,viSpk_clu1_site2);
        
        % minority of spikes has features centered other than the primary 
        % and secondary sites but ignore this for now
    %     d_mahal = mahal(mrFet1_clu1', mrFet1_clu1');
        %d = zscore(log(d));
        mrFet1_clu1 = mrFet1_clu1';
        vrDist_clu1 = madscore_(log(mahal(mrFet1_clu1, mrFet1_clu1)));
        vlExcl_clu1 = (vrDist_clu1 > thresh_mad_clu);
        if any(vlExcl_clu1)        
            S_clu.cviSpk_clu{iClu} = viSpk_clu1(~vlExcl_clu1);
            S_clu.viClu(viSpk_clu1(vlExcl_clu1)) = 0; %classify as noise cluster
            S_clu.vnSpk_clu(iClu) = numel(S_clu.cviSpk_clu{iClu});
        end
        fprintf('.');
    catch
        disperr_();
    end
end %for
fprintf('\n\ttook %0.1fs.\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [tnWav_spk, vrVrms_site] = file2spk_gt_(P, viTime_spk0)
% file loading routine. keep spike waveform (tnWav_spk) in memory
% assume that the file is chan x time format
% usage:
% [tnWav_raw, tnWav_spk, S0] = file2spk_(P)
%   
% [tnWav_raw, tnWav_spk, S0] = file2spk_(P, viTime_spk, viSite_spk)
%   construct spike waveforms from previous time markers
% 6/29/17 JJJ: Added support for the matched filter

[tnWav_spk, vnThresh_site] = deal({});    
nSamples1 = 0;  
[fid1, nBytes_file1] = fopen_(P.vcFile, 'r');
nBytes_file1 = file_trim_(fid1, nBytes_file1, P);
[nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file1, P);
t_dur1 = tic;
mnWav11_pre = [];
for iLoad1 = 1:nLoad1
    fprintf('\tProcessing %d/%d...\n', iLoad1, nLoad1);
    nSamples11 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);
    [mnWav11, vrWav_mean11] = load_file_(fid1, nSamples11, P);
    if iLoad1 < nLoad1
        mnWav11_post = load_file_preview_(fid1, P);
    else
        mnWav11_post = [];
    end    
    [viTime_spk11] = filter_spikes_(viTime_spk0, [], nSamples1 + [1, nSamples11]);
    [tnWav_spk{end+1}, vnThresh_site{end+1}] = wav2spk_gt_(mnWav11, P, viTime_spk11, mnWav11_pre, mnWav11_post);
    if iLoad1 < nLoad1, mnWav11_pre = mnWav11(end-P.nPad_filt+1:end, :); end
    nSamples1 = nSamples1 + nSamples11;
    clear mnWav11 vrWav_mean11;
end %for
fclose(fid1);   
t_dur1 = toc(t_dur1);
t_rec1 = (nBytes_file1 / bytesPerSample_(P.vcDataType) / P.nChans) / P.sRateHz;
fprintf('took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
    t_dur1, nBytes_file1/1e6, nBytes_file1/t_dur1/1e6, t_rec1/t_dur1);

% tnWav_raw = cat(3, tnWav_raw{:});
tnWav_spk = cat(3, tnWav_spk{:});
[vnThresh_site] = multifun_(@(x)cat(1, x{:}), vnThresh_site);
vrVrms_site = mean(single(vnThresh_site),1) / P.qqFactor;

end %func


%--------------------------------------------------------------------------
function [tnWav_spk, vnThresh_site] = wav2spk_gt_(mnWav1, P, viTime_spk, mnWav1_pre, mnWav1_post)
% tnWav_spk: spike waveform. nSamples x nSites x nSpikes
% trFet_spk: nSites x nSpk x nFet
% miSite_spk: nSpk x nFet
% spikes are ordered in time
% viSite_spk and viTime_spk is uint32 format, and tnWav_spk: single format
% mnWav1: raw waveform (unfiltered)
% wav2spk_(mnWav1, vrWav_mean1, P)
% wav2spk_(mnWav1, vrWav_mean1, P, viTime_spk, viSite_spk)
% 7/5/17 JJJ: accurate spike detection at the overlap region

if nargin<5, mnWav1_pre = []; end
if nargin<6, mnWav1_post = []; end

% Filter
if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)
    mnWav1 = [mnWav1_pre; mnWav1; mnWav1_post];
end
mnWav1 = int16(fft_clean_(mnWav1, P));
[mnWav2, vnWav11] = filt_car_(mnWav1, P); % filter and car

% detect spikes or use the one passed from the input (importing)
vnThresh_site = gather_(int16(mr2rms_(mnWav2, 1e5) * P.qqFactor));
nPad_pre = size(mnWav1_pre,1);
viTime_spk = viTime_spk + nPad_pre;

% reject spikes within the overlap region: problem if viTime prespecified.
if ~isempty(mnWav1_pre) || ~isempty(mnWav1_post)        
    ilim_spk = [nPad_pre+1, size(mnWav2,1) - size(mnWav1_post,1)]; %inclusive
    viKeep_spk = find(viTime_spk >= ilim_spk(1) & ilim_spk <= ilim_spk(2));
    viTime_spk = multifun_(@(x)x(viKeep_spk), viTime_spk);    
end%if
% [tnWav_spk_raw, tnWav_spk] = mn2tn_wav_(mnWav1, mnWav2, [], viTime_spk, P);
tnWav_spk = permute(gather_(mr2tr3_(mnWav2, P.spkLim, viTime_spk)), [1,3,2]);
    
% if nPad_pre > 0, viTime_spk = viTime_spk - nPad_pre; end
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_quality_(S_clu, P, viClu_update)
% 7/5/17 JJJ: Added isolation distance, L-ratio, 
% TODO: update when deleting, splitting, or merging
if nargin<3, viClu_update = []; end
t1 = tic;
fprintf('Calculating cluster quality...\n');
% [vrVmin_clu, vrVmax_clu, vrVpp_clu] = clu_amp_(S_clu.tmrWav_clu, S_clu.viSite_clu);
% [vrVmin_uV_clu, vrVmax_uV_clu, vrVpp_uV_clu] = clu_amp_(S_clu.tmrWav_raw_clu, S_clu.viSite_clu);
mrVmin_clu = squeeze_(min(S_clu.tmrWav_clu));
mrVmax_clu = squeeze_(max(S_clu.tmrWav_clu));
mrVmin_uv_clu = squeeze_(min(S_clu.tmrWav_raw_clu));
mrVmax_uv_clu = squeeze_(max(S_clu.tmrWav_raw_clu));
% tmrWav_clu = S_clu.tmrWav_clu;
% mrVmin_clu = shiftdim(min(tmrWav_clu,[],1));
% mrVmax_clu = shiftdim(max(tmrWav_clu,[],1));
vrVpp_clu = mr2vr_sub2ind_(mrVmax_clu-mrVmin_clu, S_clu.viSite_clu, 1:S_clu.nClu);
vrVmin_clu = mr2vr_sub2ind_(mrVmin_clu, S_clu.viSite_clu, 1:S_clu.nClu);
vrVpp_uv_clu = mr2vr_sub2ind_(mrVmax_uv_clu-mrVmin_uv_clu, S_clu.viSite_clu, 1:S_clu.nClu);
vrVmin_uv_clu = mr2vr_sub2ind_(mrVmin_uv_clu, S_clu.viSite_clu, 1:S_clu.nClu);
% [vrVpp_clu, ~] = max(mrVmax_clu - mrVmin_clu,[],1);
% [vrVmin_clu, viSite_min_clu] = min(mrVmin_clu,[],1);
% if ~isfield(S_clu, 'viSite_clu')
%     viSite_clu = viSite_min_clu;
% else
%     viSite_clu = S_clu.viSite_clu;
% end
% vrVpp_clu=vrVpp_clu(:);  viSite_clu=viSite_clu(:);
% [vrVpp_clu, viSite_clu, vrVmin_clu] = multifun_(@(x)x(:), vrVpp_clu, viSite_clu, vrVmin_clu);

try
    S0 = get0_();
    if isempty(S0), S0 = load0_(subsFileExt_(P.vcFile_prm, '_jrc.mat')); end
    vrVrms_site = bit2uV_(single(S0.vrThresh_site(:)) / S0.P.qqFactor, P);
%     vrSnr_clu = vrVpp_clu ./ vrVrms_site(viSite_clu);
    vrSnr_clu = abs(vrVmin_clu) ./ vrVrms_site(S_clu.viSite_clu);
    vnSite_clu = sum(bsxfun(@lt, mrVmin_clu, -vrVrms_site * S0.P.qqFactor),1)';
catch
    [vrVrms_site, vrSnr_clu, vnSite_clu] = deal([]);
    disp('no Sevt in memory.');
end
[vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = S_clu_quality2_(S_clu, P, viClu_update);

S_clu = struct_add_(S_clu, vrVpp_clu, vrSnr_clu, vrVrms_site, vnSite_clu, ...
    vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu, vrVpp_uv_clu, vrVmin_uv_clu);
fprintf('\ttook %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
function [vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = S_clu_quality2_(S_clu, P, viClu_update)
% 7/5/17 JJJ: Berenyi2013 based

global trFet_spk;
if nargin<3, viClu_update = []; end
% fprintf(2, 'S_clu_quality2_ incomplete\n');
% return; %TODO
warning off;

nSamples_2ms = round(P.sRateHz * .002);
nSamples_20ms = round(P.sRateHz * .02);
[viTime_spk, viSite_spk, viSite2_spk, cviSpk_site, cviSpk2_site] = ...
    get0_('viTime_spk', 'viSite_spk', 'viSite2_spk', 'cviSpk_site', 'cviSpk2_site');
if isempty(viClu_update)
    [vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = deal(nan(S_clu.nClu, 1));
    viClu_update = 1:S_clu.nClu;
else
    [vrIsoDist_clu, vrLRatio_clu, vrIsiRatio_clu] = get_(S_clu, 'vrIsoDist_clu', 'vrLRatio_clu', 'vrIsiRatio_clu');
    viClu_update = viClu_update(:)';
end
for iClu = viClu_update
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};    
    % Compute ISI ratio
    viTime_clu1 = viTime_spk(viSpk_clu1);
    viDTime_clu1 = diff(viTime_clu1);
    vrIsiRatio_clu(iClu) = sum(viDTime_clu1<=nSamples_2ms) ./ sum(viDTime_clu1<=nSamples_20ms);

    % Compute L-ratio an disodist (use neighboring features
    iSite_clu1 = S_clu.viSite_clu(iClu);
    % find spikes whose primary or secondary spikes reside there
    viSpk1_local = cviSpk_site{iSite_clu1}; %find(viSite_spk == iSite_clu1);
    viSpk2_local = cviSpk2_site{iSite_clu1}; %find(viSite2_spk == iSite_clu1);
    if isempty(viSpk2_local)
        mrFet12_spk = squeeze_(trFet_spk(:,1,viSpk1_local));
        viSpk12_local = viSpk1_local(:);
    else        
        mrFet12_spk = [squeeze_(trFet_spk(:,1,viSpk1_local)), squeeze_(trFet_spk(:,2,viSpk2_local))]; %squeeze_(cat(3, trFet_spk(:,1,viSpk1_local), trFet_spk(:,2,viSpk2_local)))';
        viSpk12_local = [viSpk1_local(:); viSpk2_local(:)];  
    end
    vlClu12_local = S_clu.viClu(viSpk12_local) == iClu;    
    nSpk_clu1 = sum(vlClu12_local);
    [vrLRatio_clu(iClu), vrIsoDist_clu(iClu)] = deal(nan);
    try
        vrMahal12 = mahal(mrFet12_spk', mrFet12_spk(:,vlClu12_local)');
    catch
        continue;
    end
    vrMahal12_out = vrMahal12(~vlClu12_local);
    if isempty(vrMahal12_out), continue; end
    vrLRatio_clu(iClu) = sum(1-chi2cdf(vrMahal12_out, size(mrFet12_spk,2))) / nSpk_clu1;

    % compute isolation distance
    if mean(vlClu12_local) > .5, continue; end
    sorted12 = sort(vrMahal12_out);
    vrIsoDist_clu(iClu) = sorted12(nSpk_clu1);
end %for
end %func


%--------------------------------------------------------------------------
function auto_scale_proj_time_()
% auto-scale and refgresh
S0 = get0_();

autoscale_pct = get_set_(S0.P, 'autoscale_pct', 99.5);
[hFig, S_fig] = get_fig_cache_('FigProj');
% vrY = [S_fig.hPlot0.YData(:); S_fig.hPlot1.YData(:); S_fig.hPlot2.YData(:)];
[mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, S_fig.viSites_show);
% mrAmp = [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2];
mrAmp = [mrMin1, mrMax1, mrMin2, mrMax2];
S_fig.maxAmp = quantile(mrAmp(:), autoscale_pct/100);
% S_fig.maxAmp %debug
set(hFig, 'UserData', S_fig);

[hFig, S_fig] = get_fig_cache_('FigTime');
iSite = S0.S_clu.viSite_clu(S0.iCluCopy);
[vrFet0, vrTime0] = getFet_site_(iSite);    % plot background    
[vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.iCluCopy); % plot iCluCopy
if isempty(S0.iCluPaste)
    vrFet = [vrFet0(:); vrFet1(:)];
else
    [vrFet2, vrTime2, vcYlabel, viSpk2] = getFet_site_(iSite, S0.iCluPaste); % plot iCluCopy
    vrFet = [vrFet0(:); vrFet1(:); vrFet2(:)];
end
S_fig.maxAmp = quantile(vrFet, autoscale_pct/100);
set(hFig, 'UserData', S_fig);

keyPressFcn_cell_(get_fig_cache_('FigWav'), {'j', 't'});
end %func


%--------------------------------------------------------------------------
function plot_drift_(P)
global tnWav_spk tnWav_raw trFet_spk
iShank_show = get_set_(P, 'iShank_show', 1);
vcMode_drift = get_set_(P, 'vcMode_drift', 'z'); % {'tay', 'xy', 'xya', 'x', 'y'}
vcMode_com = get_set_(P, 'vcMode_com', 'filt'); % {'fet', 'filt', 'raw', 'std'}
% Compute spike position and plot drift profile
S0 = load_cached_(P); % load cached data or from file if exists
[S_clu, viSite2_spk, vrTime_spk, viSite_spk] = get0_('S_clu', 'viSite2_spk', 'viTime_spk', 'viSite_spk');
vrTime_spk = double(vrTime_spk) / P.sRateHz;
nSites_spk = 1 + P.maxSite*2 - P.nSites_ref;
miSites_spk = single(P.miSites(:,viSite_spk));

switch lower(vcMode_com)
    case 'filt'
        %mrVp = squeeze_(single(min(tnWav_spk))) .^ 2;
%         tnWav_spk1 = meanSubt_(single(tnWav_spk),2);
        mrVp = single(squeeze_(max(tnWav_spk) - min(tnWav_spk))) .^ 2;
    case 'filtstd', mrVp = squeeze_(var(single(tnWav_spk)));
    case 'rawstd', mrVp = squeeze_(var(single(tnWav_raw)));
    case 'raw', mrVp = single(squeeze_(max(tnWav_raw) - min(tnWav_raw))) .^ 2;
    case 'fet'
        mrVp = squeeze_(trFet_spk(1:nSites_spk,1,:)) .^ 2;
        miSites_spk = miSites_spk(1:nSites_spk,:);
end


mrX_spk = reshape(P.mrSiteXY(miSites_spk,1), size(miSites_spk));
mrY_spk = reshape(P.mrSiteXY(miSites_spk,2), size(miSites_spk));
vrPosX_spk = sum(mrVp .* mrX_spk) ./ sum(mrVp);
vrPosY_spk = sum(mrVp .* mrY_spk) ./ sum(mrVp);
vrA_spk = sum(mrVp);
assignWorkspace_(vrTime_spk, vrPosX_spk, vrPosY_spk, vrA_spk);

vlSpk_shank = ismember(P.viShank_site(viSite_spk), iShank_show); %show first shank only
% vrAmp_spk = 1 ./ sqrt(single(abs(S0.vrAmp_spk)));
% vrAmp_spk = 1 ./ sqrt(sum(mrVp));
vrAmp_spk = sqrt(mean(mrVp) ./ std(mrVp)); %spatial icv
hFig_drift = create_figure_('', [0 0 .5 1], P.vcFile_prm, 1, 1);
% hFig_drift = gcf;
figure(hFig_drift); 
ax = gca(); 
hold on; 
snr_thresh_clu = get_set_(P, 'snr_thresh_clu', quantile(S_clu.vrSnr_clu, .5));
% posX_thresh = median(S_clu.vrPosX_clu);
% [vrTime_drift, vrDepth_drift] = drift_track_(S_clu, vrPosY_spk, P);

if ~isempty(S_clu)
    posX_lim = quantile(S_clu.vrPosX_clu, [.25, .75]);
    posY_lim = quantile(S_clu.vrPosY_clu, [.1, .9]);    
    for iClu = 1:S_clu.nClu
        viSpk1 = S_clu.cviSpk_clu{iClu};
        viSpk1 = viSpk1(vlSpk_shank(viSpk1));        
        vrColor1 = rand(1,3);
%         vrColor1 = 'r';
        switch vcMode_drift
            case 'tay'
                plot3(ax, vrTime_spk(viSpk1), vrAmp_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5); 
            case 'x'
                plot(ax, vrTime_spk(viSpk1), vrPosX_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5); 
            case 'y'
                plot(ax, vrTime_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5); 
            case 'xy'
                plot(ax, vrPosX_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);      
            case 'xya'
                if S_clu.vrSnr_clu(iClu) < snr_thresh_clu, continue; end
                plot3(ax, vrPosX_spk(viSpk1), vrPosY_spk(viSpk1), vrAmp_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5);                    
            case 'z'
%                 if S_clu.vrSnr_clu(iClu) < snr_thresh_clu, continue; end
%                 if S_clu.vrPosX_clu(iClu) < posX_thresh, continue; end
                posX_clu1 = S_clu.vrPosX_clu(iClu);
                posY_clu1 = S_clu.vrPosY_clu(iClu);
                if posX_clu1 < posX_lim(1) || posX_clu1 > posX_lim(2) || posY_clu1 < posY_lim(1) || posY_clu1 > posY_lim(2), continue; end
                plot(ax, vrTime_spk(viSpk1), vrPosY_spk(viSpk1), '.', 'Color', vrColor1, 'MarkerSize', 5); 
        end %switch
    end
else
    switch vcMode_drift
        case 'tay', plot3(ax, vrTime_spk, vrAmp_spk, vrPosY_spk, 'o', 'MarkerSize', 5); 
        case 'x', plot(ax, vrTime_spk, vrPosX_spk, '.', 'MarkerSize', 5); 
        case 'y', plot(ax, vrTime_spk, vrPosY_spk, '.', 'MarkerSize', 5); 
        case 'xy', plot(ax, vrPosX_spk, vrPosY_spk, '.', 'MarkerSize', 5); 
        case 'xya', plot3(ax, vrPosX_spk, vrPosY_spk, vrAmp_spk, '.', 'MarkerSize', 5); 
        case 'z'
            posX_lim = quantile(vrPosX_spk, [.25, .75]);
            posY_lim = quantile(vrPosY_spk, [.35, .65]);
            viSpk_plot = find(vrPosX_spk(:) > posX_lim(1) & vrPosX_spk(:) < posX_lim(2) & vrPosY_spk(:) > posY_lim(1) & vrPosY_spk(:) < posY_lim(2)); % & vrAmp_spk(:) < median(vrAmp_spk));
            plot(ax, vrTime_spk(viSpk_plot), vrPosY_spk(viSpk_plot), '.', 'MarkerSize', 5); 
    end %switch        
end

drawnow;
title_(ax, sprintf('Spike positions from shank %d (change using "iShank_show" paremter)', iShank_show));
xlabel(ax, 'Time (s)');
grid(ax, 'on');
axis(ax, 'tight');
switch vcMode_drift
    case 'x', ylabel(ax, 'X position (um)');
    case 'y', ylabel(ax, 'Y position (um)');
    case 'tay'
        ylabel(ax, '1/a (au)'); 
        zlabel(ax, 'z position (um)');
    case 'xya'
        xlabel(ax, 'x (um)'); 
        ylabel(ax, 'y (um)'); 
        zlabel(ax, '1/a (au)');
end %switch
mouse_figure(hFig_drift);
% linkaxes([ax_x, ax_y], 'x');
end %func


%--------------------------------------------------------------------------
% 17/9/13 JJJ: Behavior changed, if S==[], S0 is loaded
function val = get_set_(S, vcName, def_val)
% set a value if field does not exist (empty)

if isempty(S), S = get(0, 'UserData'); end
if isempty(S), val = def_val; return; end
if ~isstruct(S)
    val = []; 
    fprintf(2, 'get_set_: %s be a struct\n', inputname(1));
    return;
end
val = get_(S, vcName);
if isempty(val), val = def_val; end
end %func


%--------------------------------------------------------------------------
function S_clu = post_recenter_(S_clu, P) % recenters the features and 
% trFet_spk, viSite_spk, viSite2_spk, cviSpk_site, cviSpk2_site
global trFet_spk %nFet/site x nFet x nSpk
[viSite_spk, viSite2_spk] = get0_('viSite_spk', 'viSite2_spk');
nSites = numel(P.viSite2Chan);
% vpCentered_clu = arrayfun(@(iClu)mean(viSite_spk(S_clu.cviSpk_clu{iClu}) == S_clu.viSite_clu(iClu) | viSite2_spk(S_clu.cviSpk_clu{iClu}) == S_clu.viSite_clu(iClu)), 1:S_clu.nClu);
% find spikes not centered
for iClu = 1:S_clu.nClu
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    vlSpk_swap_clu1 = viSite2_spk(viSpk_clu1) == S_clu.viSite_clu(iClu);
    vi_swap_ = viSpk_clu1(vlSpk_swap_clu1);
    % swap indices
%     [viSite1_spk, viSite2_spk] = swap_vr_(viSite_spk, viSite2_spk, viSpk_swap_clu1);
    [viSite_spk(vi_swap_), viSite2_spk(vi_swap_)] = deal(viSite2_spk(vi_swap_), viSite_spk(vi_swap_));
    trFet_spk_clu1 = trFet_spk(:,:,vi_swap_);
    trFet_spk(:,1,vi_swap_) = trFet_spk_clu1(:,2,:);
    trFet_spk(:,2,vi_swap_) = trFet_spk_clu1(:,1,:);
end %for
cviSpk_site = arrayfun(@(iSite)find(viSite_spk == iSite), 1:nSites, 'UniformOutput', 0);
cviSpk2_site = arrayfun(@(iSite)find(viSite2_spk == iSite), 1:nSites, 'UniformOutput', 0);
S0 = set0_(viSite_spk, viSite2_spk, cviSpk_site, cviSpk2_site);
end %func


%--------------------------------------------------------------------------
function [vr1, vr2] = swap_vr_(vr1, vr2, vi_swap)
% assert: vr1 and vr2 must have the same dimensions
[vr1(vi_swap), vr2(vi_swap)] = deal(vr2(vi_swap), vr1(vi_swap));
end


%--------------------------------------------------------------------------
function [mnWav1, fGpu] = fft_clean_(mnWav, P)
fGpu = get_set_(P, 'fGpu', isGpu_(mnWav));
if ~fGpu
    mnWav1 = fft_clean(single(gather_(mnWav)), P.fft_thresh); 
else
    try
        mnWav1 = fft_clean(single(mnWav), P.fft_thresh); 
    catch
        mnWav1 = fft_clean(single(gather_(mnWav)), P.fft_thresh); 
        fGpu = 0;
    end
end
end %func


%--------------------------------------------------------------------------
function setpath_()
% Reset to the Matlab default path in case user overrided the system default functions
persistent fPathSet
if fPathSet, return; end

restoredefaultpath;
disp('Matlab path is temporarily reset to the factory default for this session.');

[vcPath_jrc, ~, ~] = fileparts(mfilename('fullpath')); % get current directory
addpath(vcPath_jrc);
fPathSet = 1;
end


%--------------------------------------------------------------------------
%Program for calculating the Adjusted Mutual Information (AMI) between
%two clusterings, tested on Matlab 7.0 (R14)
%(C) Nguyen Xuan Vinh 2008-2010
%Contact: n.x.vinh@unsw.edu.au 
%         vthesniper@yahoo.com
%--------------------------------------------------------------------------
%**Input: a contingency table T
%   OR
%        cluster label of the two clusterings in two vectors
%        eg: true_mem=[1 2 4 1 3 5]
%                 mem=[2 1 3 1 4 5]
%        Cluster labels are coded using positive integer. 
%**Output: AMI: adjusted mutual information  (AMI_max)
%
%**Note: In a prevous published version, if you observed strange AMI results, eg. AMI>>1, 
%then it's likely that in these cases the expected MI was incorrectly calculated (the EMI is the sum
%of many tiny elements, each falling out the precision range of the computer).
%However, you'll likely see that in those cases, the upper bound for the EMI will be very
%tiny, and hence the AMI -> NMI (see [3]). It is recommended setting AMI=NMI in
%these cases, which is implemented in this version.
%--------------------------------------------------------------------------
%References: 
% [1] 'A Novel Approach for Automatic Number of Clusters Detection based on Consensus Clustering', 
%       N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on 
%       Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.
% [2] 'Information Theoretic Measures for Clusterings Comparison: Is a
%	    Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J.,
%	    in Procs. the 26th International Conference on Machine Learning (ICML'09)
% [3] 'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, 
%       Normalization and Correction for Chance', N.X. Vinh, Epps, J. and
%       Bailey, J., Journal of Machine Learning Research, 11(Oct), pages
%       2837-2854, 2010

function [AMI_] = ami_(true_mem,mem)
if nargin==1
    true_mem=double(true_mem);
    T=true_mem; %contingency table pre-supplied
elseif nargin==2
    true_mem=double(true_mem);
    mem=double(mem);
    %build the contingency table from membership arrays
    R=max(true_mem);
    C=max(mem);
    n=length(mem);N=n;

    %identify & removing the missing labels
    list_t=ismember(1:R,true_mem);
    list_m=ismember(1:C,mem);
%     T=Contingency_(true_mem,mem);
    T = full(sparse(true_mem(:), mem(:), 1, max(true_mem(:)), max(mem(:))));
    T=T(list_t,list_m);
end


%-----------------------calculate Rand index and others----------
n=sum(sum(T));N=n;
C=T;
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));      %sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MIRKIN=D/t1;	    %Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;      	%Hubert 1977	%p(agree)-p(disagree)
Dri=1-RI;           %distance version of the RI
Dari=1-AR;          %distance version of the ARI
%-----------------------%calculate Rand index and others%----------

%update the true dimensions
[R C]=size(T);
if C>1 a=sum(T');else a=T';end;
if R>1 b=sum(T);else b=T;end;

%calculating the Entropies
Ha=-(a/n)*log(a/n)'; 
Hb=-(b/n)*log(b/n)';

%calculate the MI (unadjusted)
MI=0;
for i=1:R
    for j=1:C
        if T(i,j)>0, MI=MI+T(i,j)*log(T(i,j)*n/(a(i)*b(j)));end;
    end
end
MI=MI/n;

%-------------correcting for agreement by chance---------------------------
AB=a'*b;
bound=zeros(R,C);
sumPnij=0;

E3=(AB/n^2).*log(AB/n^2);

EPLNP=zeros(R,C);
LogNij=log([1:min(max(a),max(b))]/N);
for i=1:R
    for j=1:C
        sumPnij=0;
        nij=max(1,a(i)+b(j)-N);
        X=sort([nij N-a(i)-b(j)+nij]);
        if N-b(j)>X(2)
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)] [X(2)+1:N-b(j)]];
            dem=[[N-a(i)+1:N] [1:X(1)]];
        else
            nom=[[a(i)-nij+1:a(i)] [b(j)-nij+1:b(j)]];       
            dem=[[N-a(i)+1:N] [N-b(j)+1:X(2)] [1:X(1)]];
        end
%         p0=prod(nom./dem)/N;
        p0 = exp(sum(log(nom)) - sum(log(dem)))/N;
        
        sumPnij=p0;
        
        EPLNP(i,j)=nij*LogNij(nij)*p0;
        p1=p0*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);  
        
        for nij=max(1,a(i)+b(j)-N)+1:1:min(a(i), b(j))
            sumPnij=sumPnij+p1;
            EPLNP(i,j)=EPLNP(i,j)+nij*LogNij(nij)*p1;
            p1=p1*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);            
            
        end
         CC=N*(a(i)-1)*(b(j)-1)/a(i)/b(j)/(N-1)+N/a(i)/b(j);
         bound(i,j)=a(i)*b(j)/N^2*log(CC);         
    end
end

EMI_bound=sum(sum(bound));
EMI_bound_2=log(R*C/N+(N-R)*(N-C)/(N*(N-1)));
EMI=sum(sum(EPLNP-E3));

AMI_=(MI-EMI)/(max(Ha,Hb)-EMI);
NMI=MI/sqrt(Ha*Hb);


%If expected mutual information negligible, use NMI.
if abs(EMI)>EMI_bound
    fprintf('The EMI is small: EMI < %f, setting AMI=NMI',EMI_bound);
    AMI_=NMI;
end
end %func


%--------------------------------------------------------------------------
% 9/27/17 Validate parameters
function flag = validate_param_(P)
% validate P

NDIM_SORT_MAX = get_set_(P, 'nC_max', 45);

csError = {};

nSites_spk = P.maxSite * 2 + 1 - P.nSites_ref;
nFet_sort = nSites_spk * P.nPcPerChan;

if nSites_spk <= 0
    csError{end+1} = sprintf('Negative # Sites/spk. Use this formula to adjust maxSite and nSites_ref (nSites_spk = 1+maxSite*2-nSites_ref)');
end
if nFet_sort > NDIM_SORT_MAX
    csError{end+1} = sprintf('# dimensions (%d) exceeds the maximum limit for CUDA code (%d), decrease maxSite', nFet_sort, NDIM_SORT_MAX);
end

% nFet = P.
% Validate format
% if isempty(P.vcFile) && isempty(P.csFile_merge)
%     csError{end+1} = '''vcFile'' or ''csFile_merge'' must be set.';
% elseif ~isempty(P.vcFile)
%     if ~exist(P.vcFile, 'file')
%         csError{end+1} = sprintf('vcFile=''%s'' does not exist', P.vcFile);
%     end
% end
% if ~exist(P.probe_file, 'file')
%     csError{end+1} = sprintf('probe_file=''%s'' does not exist', P.probe_file);
% end
% if isempty(P.nChans), csError{end+1} = sprintf('''nChans'' must be specified.', P.probe_file); end


% Validate display


% Validate preprocess


% Validate feature


% Validate cluster


% Validate post-cluster


if isempty(csError)
    flag = true;
else
    cellfun(@(vc_)fprintf(2, '%s\n', vc_), csError);
    flag = false;
end
end %func


%--------------------------------------------------------------------------
function plot_rd_(P)
fLog_delta = 0;
max_delta = 10;
S0 = load0_(P);
S_clu = S0.S_clu;
% S_clu = get0_('S_clu');
if fLog_delta
    vrY_plot = log10(S_clu.delta);
    vcY_label = 'log10 Delta';
    vrY_plot(vrY_plot>log10(max_delta)) = log10(max_delta);
else
    vrY_plot = (S_clu.delta);
    vcY_label = 'Delta';
    vrY_plot(vrY_plot>max_delta) = max_delta;
end
figure; plot(log10(S_clu.rho), vrY_plot, '.');
xlabel('log10 Rho'); ylabel(vcY_label);
title_(P.vcFile_prm);
grid on;
end %func


%--------------------------------------------------------------------------
function vrSnr_clu = S_clu_snr_(S_clu)
global tnWav_spk
S0 = get0_();
S_clu = clu2wav_(S_clu, S0.viSite_spk, tnWav_spk);
mrVmin_clu = shiftdim(min(S_clu.tmrWav_clu,[],1));
[vrVmin_clu, viSite_clu] = min(mrVmin_clu,[],1);
vrVrms_site = single(S0.vrThresh_site(:)) / S0.P.qqFactor;
vrSnr_clu = abs(vrVmin_clu(:)) ./ bit2uV_(vrVrms_site(viSite_clu), S0.P);
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_combine_(S_clu, S_clu_redo, vlRedo_clu, vlRedo_spk)
viSpk_cluA = find(~vlRedo_spk);
[~, viCluA] = ismember(S_clu.viClu(viSpk_cluA), find(~vlRedo_clu));
S_clu.viClu(viSpk_cluA) = viCluA;
S_clu.viClu(vlRedo_spk) = S_clu_redo.viClu + max(viCluA);
end %func


%--------------------------------------------------------------------------
function delete_auto_()
% SNR based delete functionality
% Ask SNR
S0 = get0_();
[S_clu, P] = get0_('S_clu', 'P');
hFig = create_figure_('', [.5 .7 .35 .3], ['Delete Auto: ', P.vcFile]);

% Ask user which clusters to delete
plot(S_clu.vrSnr_clu(:), S_clu.vnSpk_clu(:), '.'); % show cluster SNR and spike count
xlabel('Unit SNR'); ylabel('# spikes/unit'); grid on;
set(gca,'YScale','log');
% snr_thresh = inputdlg_num_('SNR threshold: ', 'Auto-deletion based on SNR', 10); % also ask about # spikes/unit (or firing rate) @TODO
csAns = inputdlg_({'Min Unit SNR:', 'Max Unit SNR:', 'Minimum # spikes/unit'}, 'Auto-deletion based on SNR', 1, {'10', 'inf', '0'}); % also ask about # spikes/unit (or firing rate) @TODO
close(hFig);

% parse user input
if isempty(csAns), return; end
snr_min_thresh = str2double(csAns{1}); 
snr_max_thresh = str2double(csAns{2}); 
count_thresh = round(str2double(csAns{3})); 
if any(isnan([snr_min_thresh, snr_max_thresh, count_thresh]))
    msgbox_('Invalid criteria.'); return; 
end
viClu_delete = find(S_clu.vrSnr_clu(:) < snr_min_thresh | S_clu.vnSpk_clu(:) < count_thresh | S_clu.vrSnr_clu(:) > snr_max_thresh);
if isempty(viClu_delete), msgbox_('No clusters deleted.'); return; end
if numel(viClu_delete) >= S_clu.nClu, msgbox_('Cannot delete all clusters.'); return; end

% Auto delete
figure_wait_(1);
iClu_delete_post = min(S_clu.viClu) - 1;
S_clu.viClu(ismember(S_clu.viClu, viClu_delete)) = iClu_delete_post;
S_clu_new_(S_clu.viClu);
S0 = gui_update_();
figure_wait_(0);

msgbox_(sprintf('Deleted %d clusters <%0.1f SNR or <%d spikes/unit.', numel(viClu_delete), snr_min_thresh, count_thresh));
S0 = save_log_(sprintf('delete-auto <%0.1f SNR or <%d spikes/unit', snr_min_thresh, count_thresh), S0);
end %func


%--------------------------------------------------------------------------
function S = struct_select_(S, csNames, viKeep, iDimm)
if isempty(csNames), return; end
if nargin<4, iDimm = 1; end

% function test
if nargin==0, S.a=rand(10,1);S.b=rand(10,3);S.c=rand(10,3,5); csNames={'a','b','c'}; viKeep=[1,2,3,5]; iDimm=1; end

for i=1:numel(csNames)
    vcName_ = csNames{i};
    if ~isfield(S, vcName_), continue; end
    try
        val = S.(vcName_);
        if isempty(val), continue; end
        ndims_ = ndims(val);
        if ndims_==2 %find a column or row vectors
            if size(val,1)==1 || size(val,2)==1, ndims_=1; end %iscol or isrow
        end
        switch ndims_
            case 1, val = val(viKeep);
            case 2
                switch iDimm
                    case 1, val = val(viKeep,:);
                    case 2, val = val(:,viKeep);
                    otherwise, disperr_('struct_select_: invalid iDimm');
                end
            case 3
                switch iDimm
                    case 1, val = val(viKeep,:,:);
                    case 2, val = val(:,viKeep,:);
                    case 3, val = val(:,:,viKeep);
                    otherwise, disperr_('struct_select_: invalid iDimm');
                end
            otherwise, disperr_('struct_select_: invalid # of dimensions (1-3 supported)');
        end %switch
        S.(vcName_) = val;
    catch
        disperr_(sprintf('struct_select_: %s field error', vcName_));
    end
end %for
end %func


%--------------------------------------------------------------------------
function S_clu = S_clu_select_(S_clu, vlKeep_clu)
% automatically trim clusters
% 7/20/17 JJJ: auto selecting vectors and matrics
% excl vnSpk_clu, viSite_clu, vrPosX_clu, vrPosY_clu

% Quality
csNames = fieldnames(S_clu);
if isempty(csNames), return; end
viMatch_v = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_v), vlKeep_clu);

viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_t), vlKeep_clu, 3); 

viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
S_clu = struct_select_(S_clu, csNames(viMatch_c), vlKeep_clu);
end %func


%--------------------------------------------------------------------------
function flag = S_clu_valid_(S_clu)
% check the validity of the cluster by making sure the number of clusters is self-consistent.

flag = 0; %assumne invalid by default

csNames = fieldnames(S_clu);
if isempty(csNames), return; end
viMatch_v = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^v\w*_clu$'), csNames, 'UniformOutput', false));
viMatch_t = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^t\w*_clu$'), csNames, 'UniformOutput', false));
viMatch_c = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^c\w*_clu$'), csNames, 'UniformOutput', false));
viMatch_m = cellfun(@(vi)~isempty(vi), cellfun(@(cs)regexp(cs, '^m\w*_clu$'), csNames, 'UniformOutput', false));
csNames_m = csNames(viMatch_m);
csNames_m{end+1} = 'mrWavCor';
nClu = S_clu.nClu;

flag_v = all(cellfun(@(vc)numel(S_clu.(vc)) == nClu, csNames(viMatch_v))); 
flag_t = all(cellfun(@(vc)size(S_clu.(vc),3) == nClu, csNames(viMatch_t))); 
flag_c = all(cellfun(@(vc)numel(S_clu.(vc)) == nClu, csNames(viMatch_c)));
flag_m = all(cellfun(@(vc)all(size(S_clu.(vc)) == [nClu, nClu]), csNames_m));
flag = flag_v && flag_t && flag_c && flag_m;
if ~flag
    fprintf(2, 'S_clu_valid_: flag_v:%d, flag_t:%d, flag_c:%d, flag_m:%d\n', flag_v, flag_t, flag_c, flag_m);
end
end %func


%--------------------------------------------------------------------------
function [S_clu, S0] = S_clu_commit_(S_clu, vcMsg)
if nargin<2, vcMsg = ''; end
if ~S_clu_valid_(S_clu)
    fprintf(2, '%s: Cluster number is inconsistent.', vcMsg);
    S0 = get0_();
    S_clu = get_(S0, 'S_clu');
else
    S0 = set0_(S_clu);
end
end %func


%--------------------------------------------------------------------------
function [vrTime_drift, vrPos_drift] = drift_track_(S_clu, vrPosY_spk, P)
viTime_spk = double(get0_('viTime_spk'));
vrAmp_spk = get0_('vrAmp_spk');
nStep = round(get_set_(P, 'tbin_drift', 2) * P.sRateHz);
nBins = ceil(max(viTime_spk) / nStep);
vrTime_drift = ((1:nBins) * nStep - round(nStep/2)) / P.sRateHz;
viiTime_spk_drift = ceil(viTime_spk / nStep); % replaces discritize function
viClu_use = find(S_clu.vrSnr_clu(:) > quantile(S_clu.vrSnr_clu, .5)); % & S_clu.vnSpk_clu(:) > quantile(S_clu.vnSpk_clu, .5));
nClu_use = numel(viClu_use);

% spike based drift estimation
% pos_lim = quantile(vrPosY_spk, [.25, .75]);
% viSpk_use = find(vrAmp_spk < median(vrAmp_spk) & vrPosY_spk(:) > pos_lim(1) & vrPosY_spk(:) < pos_lim(2));
% vrPos_drift = accumarray(viiTime_spk_drift(viSpk_use), double(vrPosY_spk(viSpk_use)),[nBins,1],@mean,nan); %viiTime_clu1, vrPosY_spk(viSpk_clu1);

% cluster based drift estimation
mrPos_clu_drift = nan(nBins, nClu_use, 'single');
for iClu1=1:nClu_use
    iClu = viClu_use(iClu1);
    viSpk_clu1 = S_clu.cviSpk_clu{iClu};
    viiTime_clu1 = viiTime_spk_drift(viSpk_clu1);
    vrPos_clu1 = double(vrPosY_spk(viSpk_clu1));
    mrPos_clu_drift(:,iClu1) = accumarray(viiTime_clu1,vrPos_clu1,[nBins,1],@median,nan); %viiTime_clu1, vrPosY_spk(viSpk_clu1);
end
% mrPos_clu_drift = medfilt1(mrPos_clu_drift, 8);
vrActivity_clu = mean(~isnan(mrPos_clu_drift));
vrPos_drift = nanmean(mrPos_clu_drift(:,vrActivity_clu>quantile(vrActivity_clu, .9)),2);

figure; plot(vrTime_drift, vrPos_drift)
end %func


%--------------------------------------------------------------------------
function uistack_(h, vc)
try
    uistack(h, vc);
catch
end
end %func


%--------------------------------------------------------------------------
function S_fig = figWav_clu_count_(S_fig)
if nargin<1, [hFig, S_fig] = get_fig_cache_('FigWav'); end
[P, S_clu] = get0_('P', 'S_clu');
delete_(get_(S_fig, 'vhText'));
if get_(P, 'fText')
    S_fig.vhText = text_nClu_(S_clu, S_fig.hAx);
else
    S_fig.vhText = [];
end
end %func


%--------------------------------------------------------------------------
% 10/8/17 JJJ: Opens the doc file from the current JRC folder
% 7/24/17 JJJ: open the latest manual from Dropbox if the link is valid
function doc_(vcFile_doc)
% Open JRCLUST PDF documentation 
if nargin<1, vcFile_doc = 'JRCLUST manual.pdf'; end
vcFile_doc = jrcpath_(vcFile_doc);
if exist_file_(vcFile_doc)
    disp(vcFile_doc);
    open(vcFile_doc);
else
    fprintf(2, 'File does not exist: %s\n', vcFile_doc); 
end
end %func


%--------------------------------------------------------------------------
% function disperr_()
% Display error message and the error stack
function disperr_(vcMsg)
% ask user to email jrclust@vidriotech.com ? for the error ticket?
dbstack('-completenames'); % display an error stack
vcErr = lasterr();
if nargin==0
    fprintf(2, '%s\n', vcErr);
else
    fprintf(2, '%s:\n\t%s\n', vcMsg, vcErr);
end
try gpuDevice(1); disp('GPU device reset'); catch, end
end %func


%--------------------------------------------------------------------------
function P = file2struct__(vcFile_file2struct)
% James Jun 2017 May 23
% Run a text file as .m script and result saved to a struct P
% _prm and _prb can now be called .prm and .prb files

% load text file. trim and line break. remove comments.  replace 
csLines_file2struct = file2lines_(vcFile_file2struct);
csLines_file2struct = strip_comments_(csLines_file2struct);
if isempty(csLines_file2struct), P=[]; return; end

try
    eval(cell2mat(csLines_file2struct'));

    S_ws = whos(); 
    csVars = {S_ws.name};
    csVars = setdiff(csVars, {'csLines_file2struct', 'vcFile_file2struct'});
    for i=1:numel(csVars)
        eval(sprintf('a = %s;', csVars{i}));
        P.(csVars{i}) = a;
    end
catch
    disp(lasterr());
    P=[];
end
end %func



%--------------------------------------------------------------------------
% Read a text file and output cell strings separated by new lines
% 7/24/17 JJJ: Code cleanup
function csLines = file2lines_(vcFile_file2struct)
if ~exist(vcFile_file2struct, 'file')
    fprintf(2, '%s does not exist.\n', vcFile_file2struct);
    csLines = {};
    return;
end

fid = fopen(vcFile_file2struct, 'r');
csLines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

csLines = csLines{1};
end %func


%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function csLines = strip_comments_(csLines)
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)x(1)~='%', csLines));

% remove comments in the middle
for i=1:numel(csLines)
    vcLine1 = csLines{i};
    iComment = find(vcLine1=='%', 1, 'first');
    if ~isempty(iComment)
        vcLine1 = vcLine1(1:iComment-1);
    end
    vcLine1 = strrep(vcLine1, '...', '');
    if ismember(strsplit(vcLine1), {'for', 'end', 'if'})
        csLines{i} = [strtrim(vcLine1), ', ']; %add blank at the end
    else
        csLines{i} = [strtrim(vcLine1), ' ']; %add blank at the end
    end
end
% csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0);
csLines = csLines(cellfun(@(x)~isempty(x), csLines));
end %func


%--------------------------------------------------------------------------
% Strip comments from cell string
% 7/24/17 JJJ: Code cleanup
function S_bp = snr_bandpass_(P)
% compute SNR after bandpass
% read IMEC meta file and calibrate the file. Load the first bit of the file
% and bnad pass and estimate the SNR. 
% vrRmsQ_bp_site
% vrRmsQ_bp_site
% Compute SNR and apply gain correction
% show site 
% vrRmsQ_site
% vrSnr_clu
% apply a correction factor

vcFile_gain = get_set_(P, 'vcFile_gain');
nSites = numel(P.viSite2Chan);
if isempty(vcFile_gain)
    vr_Scale_uv_site = repmat(P.uV_per_bit, [nSites, 1]);
else
    try
        vrGainCorr = csvread(vcFile_gain, 1, 1);
    catch
        vrGainCorr = csvread(vcFile_gain, 0, 0);
    end
end

end %func


%--------------------------------------------------------------------------
% 7/26/17: Now accept csv file format
function vrTime_trial = loadTrial_(vcFile_trial)
% import  trial time (in seconds)
vrTime_trial = [];
try
    if ~exist(vcFile_trial, 'file'), return; end

    if matchFileExt_(vcFile_trial, '.mat')
        Strial = load(vcFile_trial);
        csFields = fieldnames(Strial);
        vrTime_trial = Strial.(csFields{1});
        if isstruct(vrTime_trial)
            vrTime_trial = vrTime_trial.times;
        end
    elseif matchFileExt_(vcFile_trial, '.csv')
        vrTime_trial = csvread(vcFile_trial);
        if isrow(vrTime_trial), vrTime_trial = vrTime_trial(:); end
    end
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
% convert a cell string to a character array separated by '\n'.
% 17/7/26 code cleanup and testing
function vc = toString_(cs)
if isempty(cs)
    vc = '';
elseif ischar(cs)
    vc = cs(:)';
elseif iscell(cs)
    vc = sprintf('%s\n', cs{:});
elseif isnumeric(cs)
    vc = sprintf('%f, ', cs(:));
end
end %func


%--------------------------------------------------------------------------
% Remove leading singular dimension
% 7/26/17 JJJ: code cleanup and testing
function val = squeeze_(val)
size_ = size(val);
if numel(size_)==2 && size_(1) == 1
    val = val';
else
    val = squeeze(val);
end
end


%--------------------------------------------------------------------------
% Return sign-preserving square root
% 7/26/17 JJJ: Code cleanup and test
function x = signsqrt_(x)
x = sign(x) .* sqrt(abs(x));
end %func


%--------------------------------------------------------------------------
% Return sign-preserving log
% 7/26/17 JJJ: Code cleanup and test
function x = signlog_(x, e)
if nargin<2, e = .0001; end
x = sign(x) .* log(abs(x + e));
end %func


%--------------------------------------------------------------------------
% Install JRCLUST 1) create user.cfg 2) compile cuda 3) compile kilosort
% 7/26/17 Code cleanup and test
function install_()
% create user.cfg
if exist('user.cfg', 'file') ~= 2
    fid = fopen('user.cfg', 'w');
    fprintf(fid, 'path_dropbox = ''C:\\Dropbox\\jrclust\\'';\n');
    fprintf(fid, 'path_backup = ''c:\\backup\\'';\n');
    fclose(fid);
    edit('user.cfg');
    msgbox_('Set path to ''path_dropbox'' and ''path_backup'' in user.cfg.');
end

compile_cuda_();
compile_ksort_();
end %func


%--------------------------------------------------------------------------
% Compile CUDA codes for JRCLUST
% 10/5/17 JJJ: Error messages converted to warning
% 7/26/17 JJJ: Code cleanup and testing
function fSuccess = compile_cuda_(S_cfg)
if nargin<1, S_cfg = read_cfg_(); end
t1 = tic;
csFiles_cu = S_cfg.csFiles_cu3; %version 3 cuda
disp('Compiling CUDA codes...');
fSuccess = 1;
S_gpu = gpuDevice(1);
if ispc()
    vcPath_nvcc = sprintf('"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v%0.1f\\bin\\nvcc"', S_gpu.ToolkitVersion);
else
    vcPath_nvcc = '/usr/local/cuda/bin/nvcc';
end

% C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\bin
for i=1:numel(csFiles_cu)
    vcFile_ = jrcpath_(csFiles_cu{i});
    vcCmd1 = sprintf('%s -ptx -m 64 -arch sm_35 "%s"', vcPath_nvcc, vcFile_);
    fprintf('\t%s\n\t', vcCmd1);
    try          
        status = system(vcCmd1);
        fSuccess = fSuccess && (status==0);        
    catch
        fprintf('\tWarning: CUDA could not be compiled: %s\n', vcFile_); 
    end
end
if ~fSuccess
    fprintf('\tWarning: CUDA could not be compiled but JRCLUST may work fine. If not, install CUDA toolkit v%0.1f and run "jrc install".\n', S_gpu.ToolkitVersion);
end
fprintf('\tFinished compiling, took %0.1fs\n', toc(t1));
end %func


%--------------------------------------------------------------------------
% Compile Kilosort code
% 10/5/17 JJJ: Error messages converted to warning
% 7/26/17 JJJ: Code cleanup and test
function fSuccess = compile_ksort_()
fSuccess = 1;
nTry = 3;
csFiles_cu = {'./kilosort/mexMPmuFEAT.cu', './kilosort/mexWtW2.cu', './kilosort/mexMPregMU.cu'};
delete ./kilosort/*.mex*;
delete ./kilosort/*.lib;
for iFile = 1:numel(csFiles_cu)
    for iTry = 1:nTry
        try
            drawnow;
            eval(sprintf('mexcuda -largeArrayDims -v %s;', csFiles_cu{iFile}));
            fprintf('Kilosort compile success for %s.\n', csFiles_cu{iFile});
            break;
        catch
            if iTry == nTry
                fprintf('\tKilosort could not be compiled: %s\n', csFiles_cu{iFile});
                fSuccess = 0;
            end
        end
    end
end
if ~fSuccess
    fprintf('\tWarning: Kilosort could not be compiled but it may work fine. If not, install Visual Studio 2013 and run "jrc install".\n');
end
end %func


%--------------------------------------------------------------------------
% Display list of toolbox and files needed
% 7/26/17 JJJ: Code cleanup and test
function [fList, pList] = disp_dependencies_()
[fList,pList] = matlab.codetools.requiredFilesAndProducts(mfilename());
if nargout==0
    disp('Required toolbox:');
    disp({pList.Name}');
    disp('Required files:');
    disp(fList');
end
end % func


%--------------------------------------------------------------------------
% @TODO: Limit the search es to .m files only
function flag = check_requirements_()
csToolbox_req = {'distrib_computing_toolbox', 'image_toolbox', 'signal_toolbox', 'statistics_toolbox'};
[fList, pList] = disp_dependencies_();

vlToolboxes = cellfun(@(vc)license('test', vc), csToolbox_req);
vlFiles = cellfun(@(vc)exist(vc, 'file')==2, fList);
flag = all([vlToolboxes(:); vlFiles(:)]); 
end %func


%--------------------------------------------------------------------------
function tr1 = fft_lowpass_(tr, fc, sRateHz)
if isempty(fc), tr1 = tr; return; end
dimm = size(tr);

mr = fft(reshape(tr, dimm(1), []));
vrFreq = fftshift((1:size(mr,1))/size(mr,1)) - .5;
vlUse = abs(vrFreq) <= fc/sRateHz/2;
mr(~vlUse,:) = 0;
tr1 = real(reshape(ifft(mr), dimm));
end %func


%--------------------------------------------------------------------------
function [cs,index] = sort_nat_(c,mode)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.4, 22 January 2011
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Set default value for mode if necessary.
if nargin < 2
	mode = 'ascend';
end

% Make sure mode is either 'ascend' or 'descend'.
modes = strcmpi(mode,{'ascend','descend'});
is_descend = modes(2);
if ~any(modes)
	error('sort_nat:sortDirection',...
		'sorting direction must be ''ascend'' or ''descend''.')
end

% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
	num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
	num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c in ascending or
% descending order, depending on mode.
[unused,index] = sortrows(comp);
if is_descend
	index = index(end:-1:1);
end
index = reshape(index,size(c));
cs = c(index);
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Expand cellstring and wild card to list of files
function csFiles = list_files_(csFiles, fSortMode)
if nargin<2
    P = get0_('P');
    fSortMode = get_set_(P, 'sort_file_merge', 1);
end
if ischar(csFiles)
    if any(csFiles=='*')
        csFiles = dir_file_(csFiles, fSortMode); %sort by dates
    else
        csFiles = {csFiles}; %make it a cell
    end
elseif iscell(csFiles)
    csFiles = cellfun(@(vc)dir_file_(vc, fSortMode), csFiles, 'UniformOutput', 0);
    csFiles = [csFiles{:}];
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Handle wild-card in makeprm case
% 7/25/17 Create a parameter file from a meta file
% Probe file can be interpreted from SpikeGLX meta file.
function [P, vcPrompt] = create_prm_file_(vcFile_bin, vcFile_prb, vcFile_template, fAsk)
if nargin<2, vcFile_prb = ''; end
if nargin<3, vcFile_template = ''; end
if nargin<4, fAsk = 1; end

[P, vcPrompt] = deal([]); 
P0 = file2struct_(jrcpath_(read_cfg_('default_prm')));  %P = defaultParam();
if ~isempty(vcFile_template)
    if exist(vcFile_template, 'file') == 2
        P0 = struct_merge_(P0, file2struct_(vcFile_template));
    else
        vcPrompt = sprintf('%s does not exist.\n', vcFile_bin);    
        fprintf(2, '%s\n', vcPrompt); 
        return;
    end    
end

if any(vcFile_bin=='*') %wild card provided
    P.csFile_merge = vcFile_bin;
    vcFile_bin = strrep(vcFile_bin, '*', '');
elseif isTextFile_(vcFile_bin)
    P.csFile_merge = vcFile_bin;
%     vcFile_bin = subsFileExt_(vcFile_bin, '.bin');
else
    if ~exist_file_(vcFile_bin)
        vcFile_bin_ = jrcpath_(vcFile_bin);
        if exist(vcFile_bin_, 'file') == 2, vcFile_bin = vcFile_bin_; end
    end
    if exist_file_(vcFile_bin)
        P.vcFile = vcFile_bin;
        P.csFile_merge = {};
    else
        vcPrompt = sprintf('%s does not exist.\n', vcFile_bin);    
        fprintf(2, '%s\n', vcPrompt); 
        return;
    end
end

% Load meta file
if isempty(P.csFile_merge)
    vcFile_meta = subsFileExt_(vcFile_bin, '.meta');
else
    csFiles_bin = filter_files_(P.csFile_merge);
    if isempty(csFiles_bin)
        vcFile_meta = '';
    else
        vcFile_meta = subsFileExt_(csFiles_bin{1}, '.meta');
    end
end
vcFile_meta = jrcpath_(vcFile_meta, 1);
P_meta = read_meta_file_(vcFile_meta);
if isempty(P_meta), P=[]; return; end

% Get the probe file if missing
if isempty(vcFile_prb)
    if isfield(P_meta.Smeta, 'imProbeOpt')
        if P_meta.Smeta.imProbeOpt > 0
            vcFile_prb = sprintf('imec3_opt%d.prb', round(P_meta.Smeta.imProbeOpt));
        end
    end
end
if isempty(vcFile_prb) % ask user
    vcFile_prb = inputdlg({'Probe file'}, 'Please specify a probe file', 1, {''});
    if isempty(vcFile_prb), P=[]; return; end
    vcFile_prb = vcFile_prb{1};
    if isempty(vcFile_prb)
        P = [];
        fprintf(2, 'You must specify a probe file (.prb)\n'); 
        return; 
    end
end

% Assign prm file name
[~,vcPostfix,~] = fileparts(vcFile_prb);
P.vcFile_prm = subsFileExt_(vcFile_bin, ['_', vcPostfix, '.prm']);
P.probe_file = vcFile_prb;
try
    S_prb = file2struct_(find_prb_(vcFile_prb));
    if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
    if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
catch
    disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
end

if exist(P.vcFile_prm, 'file') && fAsk
    vcAns = questdlg_('File already exists. Overwrite prm file?', 'Warning', 'Yes', 'No', 'No');
    if ~strcmpi(vcAns, 'Yes')
        P = [];
        vcPrompt = 'Cancelled by user.';
        return;
    end
end

% Load prb file
if isfield(P, 'template_file')
    P = struct_merge_(file2struct_(P.template_file), P);
end
P = struct_merge_(P0, P);    
P = struct_merge_(P, P_meta);    
P = struct_merge_(P, file_info_(vcFile_bin));
P.duration_file = P.nBytes_file / bytesPerSample_(P.vcDataType) / P.nChans / P.sRateHz; %assuming int16
P.version = jrc_version_();
try
    copyfile(jrcpath_(read_cfg_('default_prm')), P.vcFile_prm, 'f');
catch
    fprintf(2, 'Invalid path: %s\n', P.vcFile_prm);
    return;
end
edit_prm_file_(P, P.vcFile_prm);
vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
disp(vcPrompt);
if fAsk, edit(P.vcFile_prm); end % Show settings file
end %func


%--------------------------------------------------------------------------
% 8/4/17 JJJ: selective struct merge
% 7/31/17 JJJ: documentation and testing
function P = struct_merge_(P, P1, csNames)
% Merge second struct to first one
% P = struct_merge_(P, P_append)
% P = struct_merge_(P, P_append, var_list) : only update list of variable names

if isempty(P1), return; end
if nargin<3, csNames = fieldnames(P1); end
if ischar(csNames), csNames = {csNames}; end

for iField = 1:numel(csNames)
    vcName_ = csNames{iField};
    if isfield(P1, vcName_), P.(vcName_) = P1.(vcName_); end
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and added test ouput
function S_out = test_(vcFunc, cell_Input, nOutput, fVerbose)
% S_out = test_(vcFunc, {input1, input2, ...}, nOutput)

if nargin<2, cell_Input = {}; end
if nargin<3, nOutput = []; end
if nargin<4, fVerbose = ''; end

if isempty(nOutput), nOutput = 1; end
if ~iscell(cell_Input), cell_Input = {cell_Input}; end
if isempty(fVerbose), fVerbose = 1; end

delete_empty_files_(); % delete empty files

try
    switch nOutput
        case 0
            feval(vcFunc, cell_Input{:});
            S_out = [];
        case 1
            [out1] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1);
        case 2
            [out1, out2] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2);
        case 3
            [out1, out2, out3] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3);
        case 4
            [out1, out2, out3, out4] = feval(vcFunc, cell_Input{:});
            S_out = makeStruct_(out1, out2, out3, out4);
    end %switch
    if fVerbose
        if nOutput>=1, fprintf('[%s: out1]\n', vcFunc); disp(S_out.out1); end
        if nOutput>=2, fprintf('[%s: out2]\n', vcFunc); disp(S_out.out2); end
        if nOutput>=3, fprintf('[%s: out3]\n', vcFunc); disp(S_out.out3); end
        if nOutput>=4, fprintf('[%s: out4]\n', vcFunc); disp(S_out.out4); end
    end
catch
    disperr_();
    S_out = [];
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: documentation and testing
function P = file_info_(vcFile)
% Returns empty if file not found or multiple ones found

S_dir = dir(vcFile);
if numel(S_dir)==1
    P = struct('vcDate_file', S_dir.date, 'nBytes_file', S_dir.bytes);
else
    P = []; 
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function clear_(vcFile_prm)
% Clear JRCLUST global variables
if nargin<1, vcFile_prm = ''; end
global fDebug_ui;

% clear jrc3
clear(mfilename()); % clear persistent variables in the current file. Same as clear jrc3
clear global tnWav_spk tnWav_raw trFet_spk mnWav1 mrWav1 mnWav S_gt vcFile_prm_

set(0, 'UserData', []);   
try gpuDevice(1); catch, fprintf(2, 'GPU reset error.\n'); end
disp('Memory cleared on CPU and GPU');
fDebug_ui = [];

% clear specific parameter file. erase prm file related files
if ~isempty(vcFile_prm)    
    csFiles_del = strrep(vcFile_prm, '.prm', {'_jrc.mat', '_spkraw.jrc', '_spkwav.jrc', '_spkfet.jrc', '_log.mat'});
    delete_files_(csFiles_del);
end
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_empty_files_(vcDir)
if nargin<1, vcDir=[]; end
delete_files_(find_empty_files_(vcDir));
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and testing
function csFiles = find_empty_files_(vcDir)
% find files with 0 bytes

if nargin==0, vcDir = []; end
if isempty(vcDir), vcDir = pwd(); end
vS_dir = dir(vcDir);
viFile = find([vS_dir.bytes] == 0 & ~[vS_dir.isdir]);
csFiles = {vS_dir(viFile).name};
csFiles = cellfun(@(vc)[vcDir, filesep(), vc], csFiles, 'UniformOutput', 0);
end %func


%--------------------------------------------------------------------------
% 7/31/17 JJJ: Documentation and test
function delete_files_(csFiles, fVerbose)
% Delete list of files
% delete_files_(vcFile)
% delete_files_(csFiles)
% delete_files_(csFiles, fVerbose)

if nargin<2, fVerbose = 1; end
if ischar(csFiles), csFiles = {csFiles}; end
for iFile = 1:numel(csFiles)
    try
        if exist(csFiles{iFile}, 'file')
            delete(csFiles{iFile});
            if fVerbose
                fprintf('\tdeleted %s.\n', csFiles{iFile});
            end
        end
    catch
        disperr_();
    end
end
end %func


%--------------------------------------------------------------------------
% 9/14/17 JJJ: supports custom template.prm file
% 8/2/17 JJJ: Testing and documentation
function [vcFile_prm, vcPrompt] = makeprm_(vcFile_bin, vcFile_prb, fAsk, vcFile_template)
% Make a paramter file
% vcFile_prm = makeprm_(vcFile_bin, vcFile_prb, fAsk, vcFile_template)
% vcFile_prm = makeprm_(vcFile_bin, vcFile_template, fAsk)
global fDebug_ui
if nargin<3, fAsk = 1; end
if nargin<4, vcFile_template = ''; end

if fDebug_ui==1, fAsk = 0; end

set(0, 'UserData', []); %clear memory
if matchFileExt_(vcFile_prb, '.prm') && isempty(vcFile_template)
    vcFile_template = vcFile_prb;
    S0 = file2struct_(vcFile_template);
    vcFile_prb = S0.probe_file;    
end
[P, vcPrompt] = create_prm_file_(vcFile_bin, vcFile_prb, vcFile_template, fAsk);   
if isempty(P)
    [vcFile_prm, vcPrompt] = deal('');
    return;
end
set0_(P);
vcFile_prm = P.vcFile_prm;
end


%--------------------------------------------------------------------------
% 8/2/17 JJJ: added '.' if dir is empty
% 7/31/17 JJJ: Substitute file extension
function varargout = subsFileExt_(vcFile, varargin)
% Substitute the extension part of the file
% [out1, out2, ..] = subsFileExt_(filename, ext1, ext2, ...)

[vcDir_, vcFile_, ~] = fileparts(vcFile);
if isempty(vcDir_), vcDir_ = '.'; end
for i=1:numel(varargin)
    vcExt_ = varargin{i};    
    varargout{i} = [vcDir_, filesep(), vcFile_, vcExt_];
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function P = read_meta_file_(vcFile_meta)
% Parse meta file, ask user if meta file doesn't exist
P = [];
try
    if exist(vcFile_meta, 'file') == 2
        S_meta = read_whisper_meta_(vcFile_meta);
        P = struct('sRateHz', S_meta.sRateHz, 'uV_per_bit', S_meta.scale, 'nChans', S_meta.nChans, 'vcDataType', S_meta.vcDataType);
        %'probe_file', [S_meta.vcProbe, '.prb'], 
        P.Smeta = S_meta;    
    else
        fprintf('%s is not found. Asking users to fill out the missing info\n', vcFile_meta);
        csAns = inputdlg_({'sampling rate (Hz)', '# channels in file', 'uV/bit', 'Neuropixels option (0 if N/A)', 'Header offset (bytes)'}, 'Recording format', 1, {'30000', '385', '1','0','0'});
        if isempty(csAns), return; end
        P = struct('sRateHz', str2double(csAns{1}), 'nChans', str2double(csAns{2}), 'uV_per_bit', str2double(csAns{3}), 'imProbeOpt', str2double(csAns{4}), 'header_offset', str2double(csAns{5}));
        P.Smeta = P;
    end
catch
    disperr_('read_meta_file_');
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = read_whisper_meta_(vcFname)
% Import SpikeGLX meta file format

S = [];
viRef_imec3 = [37 76 113 152 189 228 265 304 341 380];

% Read file
if nargin < 1
    [FileName,PathName,FilterIndex] = uigetfile();
    vcFname = fullfile(PathName, FileName);
    if ~FilterIndex
        return; 
    end
end

try
    %Read Meta
    S = text2struct_(vcFname);
    S.vcDataType = 'int16'; %whisper standard
    S.ADC_bits = 16;
    S.vcProbe = '';
    
    %convert new fields to old fields    
    if isfield(S, 'niSampRate')        
        % SpikeGLX
        S.nChans = S.nSavedChans;
        S.sRateHz = S.niSampRate;
        S.rangeMax = S.niAiRangeMax;
        S.rangeMin = S.niAiRangeMin;
        S.auxGain = S.niMNGain;
        try
            S.outputFile = S.fileName;
            S.sha1 = S.fileSHA1;      
            S.vcProbe = 'imec2';
        catch
            S.outputFile = '';
            S.sha1 = [];      
        end
    elseif isfield(S, 'imSampRate')
        % IMECIII
        S.nChans = S.nSavedChans;
        S.sRateHz = S.imSampRate;
        S.rangeMax = S.imAiRangeMax;
        S.rangeMin = S.imAiRangeMin;
        S.ADC_bits = 10;  %10 bit adc but 16 bit saved
        vnIMRO = textscan(S.imroTbl, '%d', 'Delimiter', '( ),');
        vnIMRO = vnIMRO{1};
        S.auxGain = double(vnIMRO(9)); %hard code for now;
        S.auxGain_lfp = double(vnIMRO(10)); %hard code for now;
        S.vcProbe = sprintf('imec3_opt%d', vnIMRO(3));
        S.nSites = vnIMRO(4);
        S.viSites = setdiff(1:S.nSites, viRef_imec3); %sites saved
        try
            S.S_imec3 = imec3_imroTbl_(S);
        catch
            S.S_imec3 = [];
        end
    elseif isfield(S, 'sample_rate') %nick steinmetz
        S.nChans = S.n_channels_dat;
        S.sRateHz = S.sample_rate;
    end
    
     %number of bits of ADC [was 16 in Chongxi original]
    try
        S.scale = ((S.rangeMax-S.rangeMin)/(2^S.ADC_bits))/S.auxGain * 1e6;  %uVolts
    catch
        S.scale = 1;
    end    
    S.uV_per_bit = S.scale;
    if isfield(S, 'auxGain_lfp')
        S.uV_per_bit_lfp = S.scale * S.auxGain / S.auxGain_lfp;
    end
catch
    disp(lasterr);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = text2struct_(vcFname)
% Convert text file to struct

fid = fopen(vcFname, 'r');
mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
fclose(fid);
csName = mcFileMeta{1};
csValue = mcFileMeta{2};
S = struct();
for i=1:numel(csName)
    vcName1 = csName{i};
    if vcName1(1) == '~', vcName1(1) = []; end
    try         
        eval(sprintf('%s = ''%s'';', vcName1, csValue{i}));
        eval(sprintf('num = str2double(%s);', vcName1));
        if ~isnan(num)
            eval(sprintf('%s = num;', vcName1));
        end
        eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
    catch
        fprintf('%s = %s error\n', csName{i}, csValue{i});
    end
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function vcAns = questdlg_(varargin)
% Display a question dialog box
global fDebug_ui
% if get_set_([], 'fDebug_ui', 0)
if fDebug_ui == 1
    vcAns = 'Yes';
else    
    vcAns = questdlg(varargin{:});
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function varargout = get0_(varargin)
% returns get(0, 'UserData') to the workspace
% [S0, P] = get0_();
% [S0, P, S_clu] = get0_();
% [var1, var2] = get0_('var1', 'var2'); %sets [] if doesn't exist

S0 = get(0, 'UserData'); 
if ~isfield(S0, 'S_clu'), S0.S_clu = []; end
if nargin==0
    varargout{1} = S0; 
    if nargout==0, assignWorkspace_(S0); return; end
    if nargout>=1, varargout{1} = S0; end
    if nargout>=2, varargout{2} = S0.P; end
    if nargout>=3, varargout{3} = S0.S_clu; end
    return;
end
for i=1:nargin
    try                
        eval(sprintf('%s = S0.%s;', varargin{i}, varargin{i}));
        varargout{i} = S0.(varargin{i});
    catch
        varargout{i} = [];
    end
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Output message is added
% 8/2/17 JJJ: Test and documentation
function vcMsg = assignWorkspace_(varargin)
% Assign variables to the Workspace
vcMsg = {};
for i=1:numel(varargin)
    if ~isempty(varargin{i})
        assignin('base', inputname(i), varargin{i});
        vcMsg{end+1} = sprintf('assigned ''%s'' to workspace\n', inputname(i));        
    end
end
vcMsg = cell2mat(vcMsg);
if nargout==0, fprintf(vcMsg); end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function n = bytesPerSample_(vcDataType)
% Return number of bytes per data type

switch lower(vcDataType)
    case {'char', 'byte', 'int8', 'uint8'}
        n = 1;    
    case {'int16', 'uint16'}
        n = 2;
    case {'single', 'float', 'int32', 'uint32'}
        n = 4;
    case {'double', 'int64', 'uint64'}
        n = 8;
    otherwise
        n = [];
        fprintf(2, 'Unsupported data type: %s\n', vcDataType);
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function edit_prm_file_(P, vcFile_prm)
% Modify the parameter file using the variables in the P struct

csLines = file2cellstr_(vcFile_prm); %read to cell string
csLines_var = first_string_(csLines);

csName = fieldnames(P);
csValue = cellfun(@(vcField)P.(vcField), csName, 'UniformOutput',0);
for i=1:numel(csName)
    vcName = csName{i}; %find field name with 
    if isstruct(csValue{i}), continue; end %do not write struct
    vcValue = field2str_(csValue{i});
    iLine = find(strcmpi(csLines_var, vcName));
    if numel(iLine)>1 % more than one variable found
        error(['edit_prm_file_: Multiple copies of variables found: ' vcName]); 
    elseif isempty(iLine) %did not find, append
        csLines{end+1} = sprintf('%s = %s;', vcName, vcValue);
    else    
        vcComment = getCommentExpr_(csLines{iLine});        
        csLines{iLine} = sprintf('%s = %s;\t\t\t%s', vcName, vcValue, vcComment);    
    end
end
cellstr2file_(vcFile_prm, csLines);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function csLines = file2cellstr_(vcFile)
% read text file to a cell string

fid = fopen(vcFile, 'r');
csLines = {};
while ~feof(fid), csLines{end+1} = fgetl(fid); end
fclose(fid);
csLines = csLines';
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation. Added strtrim
function cs = first_string_(cs)
% Return the first string, which is typically a variable name

if ischar(cs), cs = {cs}; end

for i=1:numel(cs)
    cs{i} = strtrim(cs{i});
    if isempty(cs{i}), continue; end
    cs1 = textscan(cs{i}, '%s', 'Delimiter', {' ','='});
    cs1 = cs1{1};
    cs{i} = cs1{1};
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function vcComment = getCommentExpr_(vcExpr)
% Return the comment part of the Matlab code

iStart = strfind(vcExpr, '%');
if isempty(iStart), vcComment = ''; return; end
vcComment = vcExpr(iStart(1):end);
end %func


%---------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function vcStr = field2str_(val)
% convert a value to a strong

switch class(val)
    case {'int', 'int16', 'int32', 'uint16', 'uint32'}
        vcFormat = '%d';
    case {'double', 'float'}
        if numel(val)==1 && mod(val(1),1)==0
            vcFormat = '%d';
        else
            vcFormat = '%g';
        end
    case 'char'
        vcStr = sprintf('''%s''', val); 
        return;
    case 'cell'
        vcStr = '{';
        for i=1:numel(val)
            vcStr = [vcStr, field2str_(val{i})];
            if i<numel(val), vcStr = [vcStr, ', ']; end
        end
        vcStr = [vcStr, '}'];
        return;
    otherwise
        vcStr = '';
        fprintf(2, 'field2str_: unsupported format: %s\n', class(val));
        return;
end

if numel(val) == 1
    vcStr = sprintf(vcFormat, val);
else % Handle a matrix or array
    vcStr = '[';
    for iRow=1:size(val,1)
        for iCol=1:size(val,2)
            vcStr = [vcStr, field2str_(val(iRow, iCol))];
            if iCol < size(val,2), vcStr = [vcStr, ', ']; end
        end
        if iRow<size(val,1), vcStr = [vcStr, '; ']; end
    end
    vcStr = [vcStr, ']'];
end
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function cellstr2file_(vcFile, csLines)
% Write a cellstring to a text file

fid = fopen(vcFile, 'w');
for i=1:numel(csLines)
    fprintf(fid, '%s\n', csLines{i});
end
fclose(fid);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function S0 = set0_(varargin)
% Set(0, 'UserData')

S0 = get(0, 'UserData'); 
% set(0, 'UserData', []); %prevent memory copy operation
for i=1:nargin
    try
        S0.(inputname(i)) = varargin{i};
    catch
        disperr_();
    end
end
set(0, 'UserData', S0);
end %func


%--------------------------------------------------------------------------
% 8/2/17 JJJ: Test and documentation
function import_tsf_(vcFile_tsf)
% import tsf format (test spike file)
% create a .bin file (fTranspose = 0)
fprintf('Converting to WHISPER format (.bin and .meta)\n\t'); t1=tic;
[mnWav, Sfile] = importTSF_(vcFile_tsf);
% nChans = size(mnWav,2);
vcFile_bin = strrep(vcFile_tsf, '.tsf', '.bin');
write_bin_(vcFile_bin, mnWav);
fprintf('\n\ttook %0.1fs.\n', toc(t1));

% write .meta file
vcFile_meta = subsFileExt_(vcFile_tsf, '.meta');
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
assignWorkspace_(Sfile);
assignWorkspace_(mnWav);
fprintf('Exported to %s, %s\n', vcFile_bin, vcFile_meta);
end %func


%--------------------------------------------------------------------------
% 8/4/17 JJJ: Accepts a text file containing a list of bin files
% 7/31/17 JJJ: handles wild card combined with cells
function [csFiles_valid, viValid] = filter_files_(csFiles, fSortMode)
% Files must exist and non-zero byte
if isTextFile_(csFiles)
    csFiles = load_batch_(csFiles);
else
    if nargin<2
        P = get0_('P');
        fSortMode = get_set_(P, 'sort_file_merge', 1);
    end
    csFiles = list_files_(csFiles, fSortMode);
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
csFiles_valid = csFiles(viValid);
if ~all(vlValid)
    fprintf('Files not found:\n');
    cellfun(@(vc)fprintf('\t%s\n', vc), csFiles(~vlValid));
end
end %func


%--------------------------------------------------------------------------
% 8/4/17 JJJ: Function created, tested and documented
function flag = isTextFile_(vcFile_txt, csExt_txt)
% Return true if a file is a text file
% Accept list of text files (csExt_txt}, provided by extensions

if nargin<2, csExt_txt = {'.txt', '.batch'}; end

if ischar(vcFile_txt)
    flag = matchFileExt_(vcFile_txt, csExt_txt);
%     flag = ~isempty(regexp(lower(vcFile_txt), '.txt$'));
else
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
% 8/4/17 JJJ: Function created, tested and documented
function csLines = load_batch_(vcFile_batch)
% Load list of files in the text file, strip comments and empty lines
% Remove empty lines, comment lines starting with '%' character

csLines = file2cellstr_(vcFile_batch);
% csFiles_prm = importdata(vcFile_batch);

csLines = csLines(cellfun(@(x)~isempty(x), csLines)); % remove empty lines
csLines = cellfun(@(x)strtrim(x), csLines, 'UniformOutput', 0); % remove empty space
csLines = csLines(cellfun(@(x)~isempty(x), csLines)); %remove empty lines
csLines = csLines(cellfun(@(x)x(1)~='%', csLines)); % remove comment lines

% Removing comments that starts with "%"
% func_comment = @(vc)vc(1) == '%';
% viComment = cellfun(@(vc)func_comment(strtrim(vc)), csFiles_prm);
% csFiles_prm(viComment) = [];
end %func


%--------------------------------------------------------------------------
function [vr1, vi] = subsample_(vr, nSubsample)
vi = subsample_vr_(1:numel(vr), nSubsample);
vr1 = vr(vi);
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function [mnWav_raw, S_preview] = load_preview_(P)
% Load the subsampled dataset
% Useful for inspecting threshold and so on. filter and 
% S_preview: which file and where it came from
if ischar(P), P = loadParam_(P); end

nLoads_max_preview = get_set_(P, 'nLoads_max_preview', 30);
sec_per_load_preview = get_set_(P, 'sec_per_load_preview', 1);

% determine files to load
if isempty(P.csFile_merge)
    csFile_bin = {P.vcFile};
else
    csFile_bin = filter_files_(P.csFile_merge);
end
csFile_bin = subsample_(csFile_bin, nLoads_max_preview);

% load files
nLoads_per_file = floor(nLoads_max_preview / numel(csFile_bin));
nSamples_per_load = round(sec_per_load_preview * P.sRateHz);

% file loading loop
[mnWav_raw, cviLim_load, csFile_load] = deal({});
% [mnWav_raw, mnWav_filt] = deal({});
P.fGpu = 0;
for iFile = 1:numel(csFile_bin)
    try
        vcFile_bin_ = csFile_bin{iFile};        
        [fid_bin_, nBytes_bin_, P.header_offset] = fopen_(vcFile_bin_, 'r');
        set0_(P);
        if isempty(fid_bin_)
            fprintf(2, '.bin file does not exist: %s\n', vcFile_bin_); 
            continue; 
        end    
        fprintf('File %d/%d: %s\n\t', iFile, numel(csFile_bin), vcFile_bin_);
        nSamples_bin_ = floor(nBytes_bin_ / bytesPerSample_(P.vcDataType) / P.nChans);
        if nSamples_bin_ < nSamples_per_load % load the whole thing            
            nLoads_per_file_ = 1;
            nSamples_per_load_ = nSamples_bin_;
        else
            nLoads_per_file_ = min(nLoads_per_file, floor(nSamples_bin_ / nSamples_per_load));
            nSamples_per_load_ = nSamples_per_load;
        end
        [cvi_lim_bin, viRange_bin] = sample_skip_([1, nSamples_per_load_], nSamples_bin_, nLoads_per_file_);   
        for iLoad_ = 1:nLoads_per_file_        
            fprintf('.');
            ilim_bin_ = cvi_lim_bin{iLoad_};
            fseek_(fid_bin_, ilim_bin_(1), P);
            mnWav_raw{end+1} = load_file_(fid_bin_, diff(ilim_bin_) + 1, P);
            cviLim_load{end+1} = ilim_bin_;
            csFile_load{end+1} = vcFile_bin_;
        end   
        fprintf('\n');        
    catch
        fprintf(2, 'Loading error: %s\n', vcFile_bin_); 
    end
    fclose_(fid_bin_, 0);
end 
nLoads = numel(mnWav_raw);
mnWav_raw = cell2mat(mnWav_raw');
% if nargout>=2, mnWav_raw = cell2mat(mnWav_raw'); end
if nargout>=2
    S_preview = makeStruct_(nLoads_per_file, nLoads_max_preview, ...
        sec_per_load_preview, nSamples_per_load, nLoads, csFile_load, cviLim_load); 
end
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Tested and documented
function trWav1 = meanSubt_(trWav1, iDimm, hFunc)
% subtract mean for mr or tr
if nargin<2, iDimm = 1; end
if nargin<3, hFunc = @mean; end
if ~isa_(trWav1, 'single') && ~isa_(trWav1, 'double')
    trWav1 = single(trWav1);
end
trWav1_dim = size(trWav1);
if numel(trWav1_dim)>2, trWav1 = reshape(trWav1, trWav1_dim(1), []); end
trWav1 = bsxfun(@minus, trWav1, hFunc(trWav1, iDimm));
if numel(trWav1_dim)>2, trWav1 = reshape(trWav1, trWav1_dim); end
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Tested and documented
function flag = isa_(vr, vcClass)
% subtract mean for mr or tr
try
    if isa(vr, 'gpuArray')
        flag = strcmpi(vcClass, classUnderlying(vr));
    else
        flag = isa(vr, vcClass);
    end
catch
    disperr_();
    flag = 0;
end
end %func


%--------------------------------------------------------------------------
function xylabel_(hAx, vcXLabel, vcYLabel, vcTitle)
if isempty(hAx), hAx = gca; end
xlabel(hAx, vcXLabel);
ylabel(hAx, vcYLabel);
if nargin>=4, title_(hAx, vcTitle); end
end %func


%--------------------------------------------------------------------------
% 8/7/17 JJJ: Created, tested and documented
function vrDatenum = file_created_(csFiles, vcDateMode)
% Return the file creation time. Matlab "dir" command returns last modified time
% vcDateMode: {'lastModifiedTime', 'creationTime', 'lastAccessTime'};
if nargin<2, vcDateMode = 'creationTime'; end

if ischar(csFiles), csFiles = {csFiles}; end
vrDatenum = nan(size(csFiles));
for iFile=1:numel(csFiles)
    vcFile_ = csFiles{iFile};
    try
        vcDate = char(java.nio.file.Files.readAttributes(java.io.File(vcFile_).toPath(), vcDateMode, javaArray('java.nio.file.LinkOption', 0)));
        iStart = find(vcDate=='=', 1, 'first');
        iEnd = find(vcDate=='Z', 1, 'last');
        vcDate = vcDate(iStart+1:iEnd-1);
        vcDate(vcDate=='T') = ' '; % replace with blank
        vrDatenum(iFile) = datenum(vcDate, 'yyyy-mm-dd HH:MM:SS.FFF');
    catch
        fprintf(2, 'File not found: %d\n', vcFile_);
    end
end
% datestr(datenum1, 'yyyy-mm-dd HH:MM:SS.FFF')
end %func


%--------------------------------------------------------------------------
% 8/7/17 JJJ: Created, tested and documented
function vrDatenum = file_created_meta_(csFiles, vcDateMode)
% Read the SpikeGLX meta file for the file cration time if available
if nargin<2, vcDateMode = 'creationTime'; end
if ischar(csFiles), csFiles = {csFiles}; end
vrDatenum = nan(size(csFiles));
for iFile=1:numel(csFiles)
    vcFile_ = csFiles{iFile};
    try
        vcFile_meta_ = subsFileExt_(vcFile_, '.meta');
        if exist(vcFile_meta_, 'file')            
            S_meta_ = text2struct_(vcFile_meta_);
            vcDatenum_ = S_meta_.fileCreateTime;
            vcDatenum_(vcDatenum_=='T') = ' ';
            vrDatenum(iFile) = datenum(vcDatenum_, 'yyyy-mm-dd HH:MM:SS');
        end
    catch
        ;
    end
    if isnan(vrDatenum(iFile))
        vrDatenum(iFile) = file_created_(vcFile_, vcDateMode);
    end
end %for
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function S0 = keyPressFcn_Fig_preview_(hFig, event, S0)
if nargin<3, S0 = get0_(); end
P = get_(S0, 'P');
S_fig = get(hFig, 'UserData');
factor = 1 + 3 * key_modifier_(event, 'shift');
nSites = numel(P.viSite2Chan);

switch lower(event.Key)
    case {'uparrow', 'downarrow'}
        S_fig.maxAmp = change_amp_(event, S_fig.maxAmp, ...
            S_fig.hPlot_traces, S_fig.hPlot_traces_bad, S_fig.hPlot_traces_thresh, ...
            S_fig.hPlot_traces_spk, S_fig.hPlot_traces_spk1);
        title_(S_fig.hAx_traces, sprintf('Scale: %0.1f uV', S_fig.maxAmp));
        set(hFig, 'UserData', S_fig);     
        
    case {'leftarrow', 'rightarrow', 'home', 'end'} %navigation
        switch lower(event.Key)
            case 'leftarrow'
                nlim_bin = S_fig.nlim_bin - (S_fig.nLoad_bin) * factor; %no overlap
                if nlim_bin(1)<1
                    msgbox_('Beginning of file', 1); 
                    nlim_bin = [1, S_fig.nLoad_bin]; 
                end
            case 'rightarrow'
                nlim_bin = S_fig.nlim_bin + (S_fig.nLoad_bin + 1) * factor; %no overlap
                if nlim_bin(2) > S_fig.nSamples_bin
                    msgbox_('End of file', 1); 
                    nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
                end
            case 'home', nlim_bin = [1, S_fig.nLoad_bin]; %begging of file
            case 'end', nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin; %end of file
        end %switch
        S_fig.nlim_bin = nlim_bin;
        set(hFig, 'UserData', S_fig);
        Fig_preview_plot_(P);
        
    case 'f' %apply filter
        set(hFig, 'UserData', setfield(S_fig, 'fFilter', ~S_fig.fFilter));
        Fig_preview_plot_(P, 1);
        
     case 't' %toggle spike threshold
        set(hFig, 'UserData', setfield(S_fig, 'fThresh_spk', ~S_fig.fThresh_spk));
        Fig_preview_plot_(P, 1);
                
    case 's' %show/hide spikes
        set(hFig, 'UserData', setfield(S_fig, 'fShow_spk', ~S_fig.fShow_spk));
        Fig_preview_plot_(P, 1);          
        
    case 'g' %grid toggle on/off
        S_fig.fGrid = ~S_fig.fGrid;
        grid_([S_fig.hAx_traces, S_fig.hAx_mean, S_fig.hAx_psd, S_fig.hAx_sites], S_fig.fGrid);
        set(hFig, 'UserData', S_fig);
        menu_label_('menu_preview_view_grid', ifeq_(S_fig.fGrid, 'Hide [G]rid', 'Show [G]rid'));
        
    case 'h', msgbox_(S_fig.csHelp, 1); %help
        
    case 'r' %reset view
        P = get0_('P');
        axis(S_fig.hAx_traces, [S_fig.nlim_bin / P.sRateHz, S_fig.siteLim+[-1,1]]);
    
    case 'e' % export to workspace
        Fig_preview_export_(hFig);
        
end %switch        
end %func


%--------------------------------------------------------------------------
function Fig_preview_trange_(hFig, vc_trange, mh)
% Sets a display time range

if nargin<1, hFig = []; end
if nargin<2, vc_trange = 'custom'; end
if isempty(hFig), hFig = get_fig_cache_('Fig_preview'); end

if strcmpi(vc_trange, 'custom') % ask user input box    
    vcAns = inputdlg_('Display time range (s)', 'Time range in seconds', 1, {'.2'});
    if isempty(vcAns), return; end
    trange = str2double(vcAns{1});
else
    trange = str2double(vc_trange);
end
if isnan(trange), return; end
menu_checkbox_(mh, vc_trange);

S_fig = get(hFig, 'UserData');
P = get0_('P');
S_fig.nLoad_bin = round(trange * P.sRateHz);
nlim_bin = S_fig.nlim_bin(1) + [0, S_fig.nLoad_bin-1];
if nlim_bin(1)<1
    nlim_bin = [1, S_fig.nLoad_bin]; 
elseif nlim_bin(2) > S_fig.nSamples_bin
    nlim_bin = [-S_fig.nLoad_bin+1, 0] + S_fig.nSamples_bin;
end
S_fig.nlim_bin = nlim_bin;
set(hFig, 'UserData', S_fig);
Fig_preview_plot_(P);
end %func


%--------------------------------------------------------------------------
function Fig_preview_site_range_(hFig)
S_fig = get(hFig, 'UserData');

% Ask user and update
vcSiteFrom = num2str(S_fig.siteLim(1));
vcSiteTo = num2str(S_fig.siteLim(end));
P = get0_('P');
nSites = numel(P.viSite2Chan);
csAns = inputdlg_({'Show site from (>=1)', sprintf('Show site to (<=%d)', nSites)}, ...
    'Display site range', 1, {vcSiteFrom, vcSiteTo});
if isempty(csAns), return; end
% if isnan(site_start) || isnan(site_end), return; end;
site_start = max(str2num(csAns{1}), 1);
site_end = min(str2num(csAns{2}), nSites);
S_fig = set_(S_fig, 'siteLim', [site_start, site_end]);
set(S_fig.hAx_traces, 'YLim', S_fig.siteLim + [-1, 1]);
set(S_fig.hAx_sites, 'YLim', S_fig.siteLim + [-1, 1]);
set(hFig, 'UserData', S_fig);
end %func


%--------------------------------------------------------------------------
function menu_parent = uimenu_options_(menu_parent, csLabels, hFunc, hFig);
% create options branch in the uimenu
if nargin<4, hFig = []; end
for i=1:numel(csLabels)
    uimenu(menu_parent, 'Label', csLabels{i}, 'Callback', @(h,e)hFunc(hFig, csLabels{i}, menu_parent));
end
end %func


%--------------------------------------------------------------------------
function Fig_preview_site_thresh_(hFig)
% Site preview threshold
S_fig = get(hFig, 'UserData');
P = get0_('P');

% ask user
if isempty(S_fig.thresh_corr_bad_site)
    vc_thresh_site = '0';
else
    vc_thresh_site = num2str(S_fig.thresh_corr_bad_site);    
end
csAns = inputdlg_('Set a correlation threshold [0-1) for detecting bad sites (0 to disable)', 'thresh_corr_bad_site (set to 0 to disable)', 1, {vc_thresh_site});
if isempty(csAns), return; end
thresh_corr_bad_site = str2num(csAns{1});
if thresh_corr_bad_site>=1 || thresh_corr_bad_site<0 || isnan(thresh_corr_bad_site)
    msgbox_('Invalid range');
    return; 
end

S_fig.thresh_corr_bad_site = thresh_corr_bad_site;
set(hFig, 'UserData', S_fig);
Fig_preview_update_(hFig, S_fig, 1);

end %func


%--------------------------------------------------------------------------
function mouse_figure_(hFig, vhAx)
for iAx = 1:numel(vhAx)
    mouse_figure(hFig, vhAx(iAx));
end
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
function S_fig = preview_(P, fDebug_ui_)
% Data summary figure, interactive

global fDebug_ui
if nargin<2, fDebug_ui_ = 0; end
fDebug_ui = fDebug_ui_;
set0_(fDebug_ui);
if ischar(P), P = loadParam_(P); end
[mnWav_raw, S_preview] = load_preview_(P);
set0_(P);
nSites = size(mnWav_raw,2);

% process signal, how about common mean? 

% Bad channel metrics, do it once
mrCorr_site = corr(single(mnWav_raw));
mrCorr_site(logical(eye(size(mrCorr_site)))) = 0;
vrCorr_max_site = max(mrCorr_site);

% Create a Figure
gap = .05;
hFig = create_figure_('Fig_preview', [0 0 .5 1], P.vcFile_prm, 1, 1); %plot a summary pannel
hAx_mean = axes('Parent', hFig, 'Position',      [gap        gap         3/4-gap         1/4-gap], 'NextPlot', 'add'); 
hAx_traces = axes('Parent', hFig, 'Position',    [gap        1/4+gap     3/4-gap         3/4-gap*2], 'NextPlot', 'add'); 
hAx_sites = axes('Parent', hFig, 'Position',     [3/4+gap,   0+gap       1/4-gap*1.5     2/3-gap*2], 'NextPlot', 'add'); 
hAx_psd = axes('Parent', hFig, 'Position',       [3/4+gap,   2/3+gap     1/4-gap*1.5     1/3-gap*2], 'NextPlot', 'add'); 
linkaxes([hAx_mean, hAx_traces], 'x');

% Callback functions
Fig_preview_menu_(hFig);
set(hFig, 'KeyPressFcn', @keyPressFcn_Fig_preview_, 'BusyAction', 'cancel');
mouse_figure(hFig, hAx_traces);

% Build S_fig
[nLoads, nSamples_bin, maxAmp] = deal(S_preview.nLoads, size(mnWav_raw,1), P.maxAmp);
% nLoad_bin = S_preview.nSamples_per_load;
nLoad_bin = round(0.1 * P.sRateHz);
nlim_bin = [1, nLoad_bin];
siteLim = [1, nSites];
[vcFilter, vcCommonRef, thresh_corr_bad_site, fft_thresh, qqFactor, blank_thresh, blank_period_ms, viSiteZero] = ...
    get_(P, 'vcFilter', 'vcCommonRef', 'thresh_corr_bad_site', 'fft_thresh', 'qqFactor', 'blank_thresh', 'blank_period_ms', 'viSiteZero');
[fGrid, fFilter, fThresh_spk, fShow_spk] = deal(1, 1, 0, 1);
[vcSite_view, vcRef_view, vcPsd_view] = deal('Site correlation', 'binned', 'original');
csHelp = { ...    
    'Left/Right: change time (Shift: x4)', ...
    '[Home/End]: go to beginning/end of file', ...
    '---------', ...
    'Up/Down: change scale (Shift: x4)', ...
    'Right/Left: change time (Shift: x4)', ...
    'Zoom: Mouse wheel', ...
    '[x/y/ESC]: zoom direction', ...
    'Pan: hold down the wheel and drag', ...
    '---------', ...
    '[F]ilter toggle', ...
    'Gri[D] toggle', ...    
    };
S_fig = makeStruct_(...
    vcFilter, vcCommonRef, thresh_corr_bad_site, fft_thresh, qqFactor, blank_thresh, blank_period_ms, viSiteZero, ...
    nlim_bin, nLoad_bin, nLoads, nSamples_bin, maxAmp, ...
    mnWav_raw, vrCorr_max_site, S_preview, csHelp, ...
    hAx_mean, hAx_traces, hAx_sites, hAx_psd, ...
    fFilter, fGrid, fThresh_spk, siteLim, fShow_spk, ...
    vcSite_view, vcRef_view, vcPsd_view);

% Exit
set(hFig, 'UserData', S_fig);
drawnow;
S_fig = Fig_preview_update_(hFig, S_fig, 0);
end %func


%--------------------------------------------------------------------------
% 8/11/17 JJJ: Initial implementation
function S_fig = Fig_preview_update_(hFig, S_fig, fKeepView)
% Handle parameter update requests
% S_fig = Fig_preview_update_(hFig, S_fig, P)
% S_fig = Fig_preview_update_(hFig, S_fig, fUpdate)

if nargin<1, hFig = get_fig_cache_('Fig_preview'); end
if nargin<2, S_fig = get(hFig, 'UserData'); end
if nargin<3, fKeepView = 0; end

P = get0_('P');
nSites = numel(P.viSite2Chan);
figure_wait_(1, hFig); drawnow;
fft_thresh = S_fig.fft_thresh;
S_fig.mnWav_clean = int16(fft_clean_(S_fig.mnWav_raw, struct_add_(P, fft_thresh))); % fft filter

% Find bad sites
if S_fig.thresh_corr_bad_site > 0
    S_fig.vlSite_bad = S_fig.vrCorr_max_site < S_fig.thresh_corr_bad_site;
    S_fig.viSite_bad = find(S_fig.vlSite_bad);
elseif ~isempty(S_fig.viSiteZero)
    S_fig.viSite_bad = S_fig.viSiteZero;    
    S_fig.vlSite_bad = false(size(S_fig.vrCorr_max_site));
    S_fig.vlSite_bad(S_fig.viSiteZero) = 1;
    S_fig.viSiteZero = [];
else
    S_fig.vlSite_bad = false(size(S_fig.vrCorr_max_site));
    S_fig.viSite_bad = [];
end

% Perform filter fft_thresh
P_ = set_(P, 'vcCommonRef', 'none', 'fGpu', 0, 'vcFilter', S_fig.vcFilter, ...
    'blank_period_ms', S_fig.blank_period_ms, 'blank_thresh', S_fig.blank_thresh, 'fParfor', 0);
mnWav_filt = filt_car_(S_fig.mnWav_clean, P_);
% if strcmpi(S_fig.vcCommonRef, 'median')
%     vrWav_filt_mean = median(mnWav_filt(:,~S_fig.vlSite_bad), 2);
% else
%     vrWav_filt_mean = mean(mnWav_filt(:,~S_fig.vlSite_bad), 2);
% end
vrWav_filt_mean = mr2ref_(mnWav_filt, S_fig.vcCommonRef, S_fig.viSite_bad);
if ~strcmpi(S_fig.vcCommonRef, 'none')
    mnWav_filt = bsxfun(@minus, mnWav_filt, int16(vrWav_filt_mean));
end
S_fig.vrWav_filt_mean = madscore_(mean(vrWav_filt_mean, 2)); % Save in MAD unit
S_fig.mnWav_filt = mnWav_filt;
% reference threshold. determine 
% vrWav_filt_mean = mean(mnWav_filt, 2) * P.uV_per_bit; % @TODO: save in MAD unit

% vrPower_psd = abs(mean(fft(S_fig.mnWav_raw(:,~S_fig.vlSite_bad)), 2));
[mrPower_psd, S_fig.vrFreq_psd] = psd_(S_fig.mnWav_raw(:,~S_fig.vlSite_bad), P.sRateHz, 4);
[mrPower_clean_psd] = psd_(S_fig.mnWav_clean(:,~S_fig.vlSite_bad), P.sRateHz, 4);
[S_fig.vrPower_psd, S_fig.vrPower_clean_psd] = multifun_(@(x)mean(x,2), mrPower_psd, mrPower_clean_psd);

% Apply threshold and perform spike detection
vrRmsQ_site = mr2rms_(mnWav_filt, 1e5);
vnThresh_site = int16(vrRmsQ_site * S_fig.qqFactor);
vnThresh_site(S_fig.vlSite_bad) = nan; % shows up as 0 for int16
S_fig.mlWav_thresh = bsxfun(@lt, mnWav_filt, -abs(vnThresh_site)); %negative threshold crossing
S_fig.mlWav_thresh(:, S_fig.vlSite_bad) = 0;
S_fig.vnThresh_site = vnThresh_site;

% Spike detection
% P_.fMerge_spk = 0;
[vlKeep_ref, S_fig.vrMad_ref] = car_reject_(vrWav_filt_mean, P_);
[S_fig.viTime_spk, S_fig.vnAmp_spk, viSite_spk] = detect_spikes_(mnWav_filt, vnThresh_site, vlKeep_ref, P_);
t_dur = size(mnWav_filt,1) / P.sRateHz;
S_fig.vrEventRate_site = hist(viSite_spk, 1:nSites) / t_dur; % event count
S_fig.vrEventSnr_site = abs(single(arrayfun(@(i)median(S_fig.vnAmp_spk(viSite_spk==i)), 1:nSites))) ./ vrRmsQ_site;
S_fig.vlKeep_ref = vlKeep_ref;
S_fig.viSite_spk = viSite_spk;
% Spike stats: such as # sites/event over threshold

% Exit
set(hFig, 'UserData', S_fig);
[hFig, S_fig] = Fig_preview_plot_(P, fKeepView);
end %func


%--------------------------------------------------------------------------
function [mrPower, vrFreq] = psd_(mr, Fs, nSkip)
% power spectrum
if nargin<3, nSkip = 1; end
if nSkip>1
    mr = mr(1:nSkip:end,:);
%     Fs = Fs / nSkip;
end
n = size(mr,1);
n1 = round(n/2);
mrPower = fft(meanSubt_(single(mr)));
mrPower = pow2db_(abs(mrPower(2:n1+1, :))) / n;

if nargout>=2, vrFreq = Fs*(1:n1)'/n; end
end %func


%--------------------------------------------------------------------------
% 8/6/17 JJJ: Initial implementation, documented and tested
% @TODO: run spike detection and show detected spikes, or below threshold values
function [hFig, S_fig] = Fig_preview_plot_(P, fKeepView)

if nargin<1, P = []; end
if nargin<2, fKeepView = 0; end
hWait = msgbox_('Plotting...', 0, 1);
if isempty(P), P = get0_('P'); end
[hFig, S_fig] = get_fig_cache_('Fig_preview');

figure_wait_(1, hFig);
nSites = size(S_fig.mnWav_filt,2);
viPlot = S_fig.nlim_bin(1):S_fig.nlim_bin(2);
vrTime_sec = viPlot / P.sRateHz;
tlim_sec = (S_fig.nlim_bin + [-1 1]) / P.sRateHz;


%-----
% Mean plot
if isempty(get_(S_fig, 'hPlot_mean'))
    S_fig.hPlot_mean = plot(S_fig.hAx_mean, nan, nan, 'k');
    S_fig.hPlot_mean_thresh = plot(S_fig.hAx_mean, nan, nan, 'r');
end
switch S_fig.vcRef_view
    case 'original'
        vrY_ = S_fig.vrWav_filt_mean(viPlot);        
    case 'binned'
        vrY_ = S_fig.vrMad_ref(viPlot);
    otherwise, disperr_();
end
set(S_fig.hPlot_mean, 'XData', vrTime_sec, 'YData', vrY_);
set(S_fig.hAx_mean, 'YLim', [-50, 50]);
xylabel_(S_fig.hAx_mean, 'Time (sec)', sprintf('Common ref. (MAD, %s)', S_fig.vcRef_view));
fThresh_ref = strcmpi(S_fig.vcRef_view, 'binned') && ~isempty(S_fig.blank_thresh) && S_fig.blank_thresh ~= 0;
if fThresh_ref
    set(S_fig.hPlot_mean_thresh, 'XData', vrTime_sec([1,end]), 'YData', repmat(S_fig.blank_thresh, [1,2]));
else
    hide_plot_(S_fig.hPlot_mean_thresh);
end

%-----
% Traces plot
if isempty(get_(S_fig, 'hPlot_traces'))
    S_fig.hPlot_traces = plot(S_fig.hAx_traces, nan, nan, 'Color', [1 1 1]*.5);
    S_fig.hPlot_traces_spk = plot(S_fig.hAx_traces, nan, nan, 'm.-', 'LineWidth', 1.5);
    S_fig.hPlot_traces_spk1 = plot(S_fig.hAx_traces, nan, nan, 'ro');
    S_fig.hPlot_traces_thresh = plot(S_fig.hAx_traces, nan, nan, 'm:');
    S_fig.hPlot_traces_bad = plot(S_fig.hAx_traces, nan, nan, 'r');
end
if S_fig.fFilter    
    vcFilter = S_fig.vcFilter;
    mrWav_ = bit2uV_(S_fig.mnWav_filt(viPlot,:), struct_add_(P, vcFilter));
else
    mrWav_ = meanSubt_(single(S_fig.mnWav_clean(viPlot,:)) * P.uV_per_bit);
end
multiplot(S_fig.hPlot_traces, S_fig.maxAmp, vrTime_sec, mrWav_, 1:nSites);
if ~isempty(S_fig.viSite_bad)
    multiplot(S_fig.hPlot_traces_bad, S_fig.maxAmp, vrTime_sec, mrWav_(:,S_fig.viSite_bad), S_fig.viSite_bad);
else
    hide_plot_(S_fig.hPlot_traces_bad);
end
if S_fig.fShow_spk
    vlSpk_ = S_fig.viTime_spk >= S_fig.nlim_bin(1) & S_fig.viTime_spk <= S_fig.nlim_bin(end);
    viTime_spk_ = single(S_fig.viTime_spk(vlSpk_)-S_fig.nlim_bin(1)+1);
    vrTime_spk_ = single(S_fig.viTime_spk(vlSpk_)) / P.sRateHz;
    viSite_spk_ = single(S_fig.viSite_spk(vlSpk_));
else
    viTime_spk_ = [];
end
if isempty(viTime_spk_)
    hide_plot_(S_fig.hPlot_traces_spk1);
    menu_label_('menu_preview_view_spike', 'Show [S]pikes');
else
    multiplot(S_fig.hPlot_traces_spk1, S_fig.maxAmp, vrTime_spk_, mr2vr_sub2ind_(mrWav_, viTime_spk_, viSite_spk_), viSite_spk_, 1);
    menu_label_('menu_preview_view_spike', 'Hide [S]pikes');
end
vrThresh_site_uV = bit2uV_(-S_fig.vnThresh_site(:), setfield(P, 'vcFilter', S_fig.vcFilter));
vrThresh_site_uV(S_fig.viSite_bad) = nan;
if S_fig.fThresh_spk && S_fig.fFilter
    multiplot(S_fig.hPlot_traces_thresh, S_fig.maxAmp, vrTime_sec([1,end,end])', repmat(vrThresh_site_uV, [1,3])');
    multiplot(S_fig.hPlot_traces_spk, S_fig.maxAmp, vrTime_sec, ...
        mr_set_(mrWav_, ~S_fig.mlWav_thresh(viPlot,:), nan)); %show spikes
    menu_label_('menu_preview_view_threshold', 'Hide spike [T]threshold');
else
    hide_plot_([S_fig.hPlot_traces_thresh, S_fig.hPlot_traces_spk]);
    menu_label_('menu_preview_view_threshold', 'Show spike [T]hreshold');
end
xylabel_(S_fig.hAx_traces, '', 'Site #');
vcFilter_ = ifeq_(S_fig.fFilter, sprintf('Filter=%s', S_fig.vcFilter), 'Filter off');
set(hFig, 'Name', sprintf('%s; %s; CommonRef=%s', P.vcFile_prm, vcFilter_, S_fig.vcCommonRef));

title_(S_fig.hAx_traces, sprintf('Scale: %0.1f uV', S_fig.maxAmp));
menu_label_('menu_preview_view_filter', ifeq_(S_fig.fFilter, 'Show raw traces [F]', 'Show [F]iltered traces'));
menu_label_('menu_preview_view_grid', ifeq_(S_fig.fGrid, 'Hide [G]rid', 'Show [G]rid'));
if ~fKeepView
    set(S_fig.hAx_traces, 'YTick', 1:nSites, 'YLim', S_fig.siteLim + [-1,1], 'XLim', tlim_sec);
end


%-----
% Site plot
if isempty(get_(S_fig, 'hPlot_site'))
    S_fig.hPlot_site = barh(S_fig.hAx_sites, nan, nan, 1); 
    S_fig.hPlot_site_bad = barh(S_fig.hAx_sites, nan, nan, 1, 'r'); 
    S_fig.hPlot_site_thresh = plot(S_fig.hAx_sites, nan, nan, 'r'); 
end
switch S_fig.vcSite_view
    case 'Site correlation', vrPlot_site = S_fig.vrCorr_max_site;
    case 'Spike threshold', vrPlot_site = single(S_fig.vnThresh_site);
    case 'Event rate (Hz)', vrPlot_site = S_fig.vrEventRate_site;
    case 'Event SNR (median)', vrPlot_site = S_fig.vrEventSnr_site;
end
set(S_fig.hPlot_site, 'XData', 1:nSites, 'YData', vrPlot_site); %switch statement
xylabel_(S_fig.hAx_sites, S_fig.vcSite_view, 'Site #');
set(S_fig.hAx_sites, 'YLim', S_fig.siteLim + [-1,1]);
if isempty(S_fig.thresh_corr_bad_site) || ~strcmpi(S_fig.vcSite_view, 'Site correlation')
    hide_plot_(S_fig.hPlot_site_thresh);
else
    set(S_fig.hPlot_site_thresh, 'XData', S_fig.thresh_corr_bad_site *[1,1], 'YData', [0, nSites+1]);
end
if ~isempty(S_fig.viSite_bad);
    vrPlot_site_bad = vrPlot_site;
    vrPlot_site_bad(~S_fig.vlSite_bad) = 0;    
    set(S_fig.hPlot_site_bad, 'XData', 1:nSites, 'YData', vrPlot_site_bad); %switch statement
else
    hide_plot_(S_fig.hPlot_site_bad);
end
title_(S_fig.hAx_sites, sprintf('thresh_corr_bad_site=%0.4f', S_fig.thresh_corr_bad_site));


%-----
% PSD plot
if isempty(get_(S_fig, 'hPlot_psd'))
    S_fig.hPlot_psd = plot(S_fig.hAx_psd, nan, nan, 'k');
    S_fig.hPlot_clean_psd = plot(S_fig.hAx_psd, nan, nan, 'g');
    S_fig.hPlot_psd_thresh = plot(S_fig.hAx_psd, nan, nan, 'r');
end
set(S_fig.hPlot_psd, 'XData', S_fig.vrFreq_psd, 'YData', S_fig.vrPower_psd);
set(S_fig.hPlot_clean_psd, 'XData', S_fig.vrFreq_psd, 'YData', S_fig.vrPower_clean_psd);
xylabel_(S_fig.hAx_psd, 'Frequency (Hz)', 'Power [dB]', 'TODO: before and after cleaning');
set(S_fig.hAx_psd, 'XLim', [0, P.sRateHz/2]);
title_(S_fig.hAx_psd, sprintf('fft_thresh=%s', num2str(S_fig.fft_thresh)));

grid_([S_fig.hAx_traces, S_fig.hAx_mean, S_fig.hAx_sites, S_fig.hAx_psd], S_fig.fGrid);

% Exit
set(hFig, 'UserData', S_fig);
figure_wait_(0, hFig);
close_(hWait);
end %func


%--------------------------------------------------------------------------
function mr = mr_set_(mr, ml, val)
mr(ml)=val;
end %func


%--------------------------------------------------------------------------
function hide_plot_(vhPlot)
for i=1:numel(vhPlot), set(vhPlot(i), 'XData', nan, 'YData', nan); end
end %func


%--------------------------------------------------------------------------
function Fig_preview_ref_(hFig, vcMode, hMenu)
% Set reference types
S_fig = get(hFig, 'UserData');
S_fig.vcCommonRef = vcMode;
S_fig = Fig_preview_update_(hFig, S_fig, 1);
menu_checkbox_(hMenu, vcMode);
end %func


%--------------------------------------------------------------------------
function Fig_preview_filter_(hFig, vcMode, hMenu)
% Set reference types
S_fig = get(hFig, 'UserData');
menu_checkbox_(hMenu, vcMode);
S_fig = Fig_preview_update_(hFig, set_(S_fig, 'vcFilter', vcMode, 'fFilter', 1), 1);
end %func


%--------------------------------------------------------------------------
function Fig_preview_fft_thresh_(hFig)
S_fig = get(hFig, 'UserData');
P = get0_('P');

% ask user
if isempty(S_fig.fft_thresh)
    vc_ = '0';
else
    vc_ = num2str(S_fig.fft_thresh);    
end
vcAns = inputdlg_('Set a threshold for FFT cleanup (0 to disable)', 'fft_thresh (20 recommended)', 1, {vc_});
if isempty(vcAns), return; end
fft_thresh = str2double(vcAns{1});
if isnan(fft_thresh), return; end

S_fig.fft_thresh = fft_thresh;
S_fig = Fig_preview_update_(hFig, S_fig, 1);
end %func


%--------------------------------------------------------------------------
function grid_(vhAx, fGrid)
if isempty(vhAx), vhAx = gca; end
for iAx = 1:numel(vhAx)
    if ischar(fGrid)
        grid(vhAx(iAx), fGrid);
    else
        grid(vhAx(iAx), ifeq_(fGrid, 'on', 'off'));
    end
end
end %func


%--------------------------------------------------------------------------
function Fig_preview_site_plot_(hFig, vcMode, hMenu)
S_fig = get(hFig, 'UserData');
S_fig.vcSite_view = vcMode;
menu_checkbox_(hMenu, vcMode);
set(hFig, 'UserData', S_fig);
Fig_preview_plot_([], 1);
end %func


%--------------------------------------------------------------------------
% 8/12/17 JJJ: Created
function Fig_preview_psd_plot_(hFig, vcMode, hMenu)
S_fig = get(hFig, 'UserData');
switch vcMode
%     case 'Power'
%     case 'Detrended'
    case 'Linear', set(S_fig.hAx_psd, 'XScale', 'linear');
    case 'Log', set(S_fig.hAx_psd, 'XScale', 'log');
    otherwise, disperr_(vcMode);
end %switch
menu_checkbox_(hMenu, vcMode);
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function Fig_preview_export_(hFig)
% Export raw and filtered traces

S_fig = get(hFig, 'UserData');
[mnWav_filt, mnWav_raw, S_preview] = deal(S_fig.mnWav_filt, S_fig.mnWav_raw, S_fig.S_preview);
assignWorkspace_(mnWav_filt, mnWav_raw, S_preview);
msgbox_('Exported to workspace: Raw traces (mnWav_raw), filtered traces (mnWav_filt), Preview info (S_preview)');
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function Fig_preview_save_prm_(hFig)
% Update to a parameter file, show a preview window
% List of parameters to update. Set threshold: save to a known location and load 
P = get0_('P');
S_fig = get(hFig, 'UserData');

% Select variables to export
[fft_thresh, qqFactor, vlSite_bad, vcFilter, vcCommonRef, blank_thresh, blank_period_ms] = ...
    get_(S_fig, 'fft_thresh', 'qqFactor', 'vlSite_bad', 'vcFilter', 'vcCommonRef', 'blank_thresh', 'blank_period_ms');
viSiteZero = find(vlSite_bad);
P_update = makeStruct_(fft_thresh, viSiteZero, qqFactor, vcFilter, vcCommonRef, blank_thresh, blank_period_ms);

% Preview variables in the edit box
vcUpdate = struct2str_(P_update);
csAns = inputdlg_(P.vcFile_prm, 'Update confirmation', 16, {vcUpdate}, struct('Resize', 'on'));
if isempty(csAns), return; end
[P_update, vcErr] = str2struct_(csAns{1});
if isempty(P_update)
    msgbox_(vcErr);
    return; 
end
P = get0_('P');
edit_prm_file_(P_update, P.vcFile_prm);
% P = struct_merge_(P, P_update);
% set0_(P);
edit(P.vcFile_prm);
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function vc = struct2str_(S)
if ~isstruct(S), vc = ''; return; end

csVars = fieldnames(S);
vc = '';
for iVar = 1:numel(csVars)
    vc = sprintf('%s%s = %s;', vc, csVars{iVar}, field2str_(S.(csVars{iVar})));
    if iVar<numel(csVars), vc = sprintf('%s\n', vc); end
end %for
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function [S, vcErr] = str2struct_(vc)
% Numerical expressions only for now

vc = vc';
vc = vc(:)';
vcErr = '';
cc_ = textscan(strtrim(vc), '%s%s', 'Delimiter', '=;',  'ReturnOnError', false);
csName = cc_{1};
csValue = cc_{2};

S = struct();
for i=1:numel(csName)
    vcName_ = strtrim(csName{i});
    if vcName_(1) == '~', vcName_(1) = []; end
    try         
        eval(sprintf('%s = %s;', vcName_, csValue{i}));
        eval(sprintf('num = str2double(%s);', vcName_));
        if ~isnan(num)
            eval(sprintf('%s = num;', vcName_));
        end
        eval(sprintf('S = setfield(S, ''%s'', %s);', vcName_, vcName_));
    catch
        vcErr = lasterr();
        S = []; 
        return;
%         fprintf('%s = %s error\n', csName{i}, csValue{i});
    end
end

end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function hMenu = menu_label_(vcTag, vcLabel)
hMenu = [];
try
    if ischar(vcTag)
        hMenu = findobj('Tag', vcTag, 'Type', 'uimenu');
    else
        hMenu = vcTag;
    end
    set(hMenu, 'Label', vcLabel); %figure property
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function hMenu = menu_checkbox_(vcTag, vcLabel, fUncheckRest)
% Check the label selected and uncheck the rest
if nargin<3, fUncheckRest = 1; end
hMenu = [];
try
    if ischar(vcTag)
        hMenu = findobj('Tag', vcTag, 'Type', 'uimenu');
    else
        hMenu = vcTag;
    end
    % Find children
    vhMenu_child = hMenu.Children;
    for iMenu = 1:numel(vhMenu_child)
        hMenu_child_ = vhMenu_child(iMenu);
        vcLabel_ = get(hMenu_child_, 'Label');
        if strcmpi(vcLabel_, vcLabel)
            set(hMenu_child_, 'Checked', 'on');
        elseif fUncheckRest
            set(hMenu_child_, 'Checked', 'off');
        end
    end %for
catch
    disperr_();
end
end %func


%--------------------------------------------------------------------------
function Fig_preview_spk_thresh_(hFig)
% update spike threshold qqFactor
S_fig = get(hFig, 'UserData');
P = get0_('P');

% ask user
vcAns = inputdlg_('Set a spike detection threshold (qqFactor)', 'Spike detection threshold', 1, {num2str(S_fig.qqFactor)});
if isempty(vcAns), return; end
qqFactor = str2num(vcAns{1});
if isnan(qqFactor) || isempty(qqFactor), return; end

S_fig = Fig_preview_update_(hFig, setfield(S_fig, 'qqFactor', qqFactor), 1);
end %func


%--------------------------------------------------------------------------
% 8/16/17 JJJ: Updated
function Fig_preview_menu_(hFig)
P = get0_('P');

set(hFig, 'MenuBar','None');
mh_file = uimenu(hFig, 'Label','File'); 
uimenu(mh_file, 'Label', sprintf('Save to %s', P.vcFile_prm), 'Callback', @(h,e)Fig_preview_save_prm_(hFig));
uimenu(mh_file, 'Label', '[E]xport to workspace', 'Callback', @(h,e)Fig_preview_export_(hFig));
uimenu(mh_file, 'Label', 'Save spike detection threshold', 'Callback', @(h,e)Fig_preview_save_threshold_(hFig));

% Edit menu
mh_edit = uimenu(hFig, 'Label','Edit'); 

uimenu(mh_edit, 'Label', 'Bad site threshold', 'Callback', @(h,e)Fig_preview_site_thresh_(hFig));
uimenu(mh_edit, 'Label', 'Spike detection threshold', 'Callback', @(h,e)Fig_preview_spk_thresh_(hFig));

mh_edit_filter = uimenu(mh_edit, 'Label', 'Filter mode');
uimenu_options_(mh_edit_filter, {'sgdiff', 'bandpass'}, @Fig_preview_filter_, hFig);
menu_checkbox_(mh_edit_filter, get_filter_(P));

mh_edit_ref = uimenu(mh_edit, 'Label', 'Reference mode');
uimenu_options_(mh_edit_ref, {'none', 'mean', 'median'}, @Fig_preview_ref_, hFig); % @TODO: local mean
menu_checkbox_(mh_edit_ref, P.vcCommonRef);

uimenu(mh_edit, 'Label', 'Common reference threshold', 'Callback', @(h,e)Fig_preview_ref_thresh_(hFig));

uimenu(mh_edit, 'Label', 'FFT cleanup threshold', 'Callback', @(h,e)Fig_preview_fft_thresh_(hFig));


% View menu
mh_view = uimenu(hFig, 'Label','View'); 

mh_view_trange = uimenu(mh_view, 'Label', 'Display time range (s)');
uimenu_options_(mh_view_trange, {'0.05', '0.1', '0.2', '0.5', '1', '2', '5', 'Custom'}, @Fig_preview_trange_, hFig);
menu_checkbox_(mh_view_trange, '0.1');
uimenu(mh_view, 'Label', 'Display site range', 'Callback', @(h,e)Fig_preview_site_range_(hFig));

uimenu(mh_view, 'Label', 'Show raw traces [F]', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'f'), 'Tag', 'menu_preview_view_filter');
uimenu(mh_view, 'Label', 'Show spike [T]hreshold', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 't'), 'Tag', 'menu_preview_view_threshold');
uimenu(mh_view, 'Label', 'Hide [S]pikes', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 's'), 'Tag', 'menu_preview_view_spike');
uimenu(mh_view, 'Label', 'Hide [G]rid', 'Callback', @(h,e)keyPressFcn_cell_(hFig, 'g'), 'Tag', 'menu_preview_view_grid');

mh_view_site = uimenu(mh_view, 'Label', 'Site view');
uimenu_options_(mh_view_site, {'Site correlation', 'Spike threshold', 'Event rate (Hz)', 'Event SNR (median)'}, @Fig_preview_site_plot_, hFig);
menu_checkbox_(mh_view_site, 'Site correlation');

mh_view_ref = uimenu(mh_view, 'Label', 'Reference view');
uimenu_options_(mh_view_ref, {'original', 'binned'}, @Fig_preview_view_ref_, hFig);
menu_checkbox_(mh_view_ref, 'original');

mh_view_freq = uimenu(mh_view, 'Label', 'Frequency scale');
uimenu_options_(mh_view_freq, {'Linear', 'Log'}, @Fig_preview_psd_plot_, hFig);
menu_checkbox_(mh_view_freq, 'Linear');

% mh_view_psd = uimenu(mh_view, 'Label', 'PSD view');
% uimenu_options_(mh_view_psd, {'original', 'detrended'}, @Fig_preview_view_psd_, hFig);
% menu_checkbox_(mh_view_psd, 'original');

% 'Power', 'Detrended', 
end %func


%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_save_threshold_(hFig)
% export S_fig.vnThresh and sets vcFile_thresh
P = get0_('P');
S_fig = get(hFig, 'UserData');
vnThresh_site = S_fig.vnThresh_site;
vcFile_thresh = strrep(P.vcFile_prm, '.prm', '_thresh.mat');
save(vcFile_thresh, 'vnThresh_site'); % also need to store filter values?
P.vcFile_thresh = vcFile_thresh;
set0_(P);
edit_prm_file_(P, P.vcFile_prm);

msgbox_(sprintf('Saved to %s and updated %s (vcFile_thresh)', ...
    vcFile_thresh, P.vcFile_prm));
end %func


%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_view_psd_(hFig, vcMode, hMenu)
% export S_fig.vnThresh and sets vcFile_thresh
S_fig = get(hFig, 'UserData');
S_fig.vcPsd_view = vcMode;
menu_checkbox_(hMenu, vcMode);
set(hFig, 'UserData', S_fig);
Fig_preview_plot_([], 1);
end %func


%--------------------------------------------------------------------------
% 8/16/17 JJJ: created
function Fig_preview_view_ref_(hFig, vcMode, hMenu)
S_fig = get(hFig, 'UserData');
S_fig.vcRef_view = vcMode;
menu_checkbox_(hMenu, vcMode);
set(hFig, 'UserData', S_fig);
Fig_preview_plot_([], 1);
end %func


%--------------------------------------------------------------------------
% 8/17/17 JJJ: created
function Fig_preview_ref_thresh_(hFig)
S_fig = get(hFig, 'UserData');

% Ask user and update
vc_blank_thresh = num2str(S_fig.blank_thresh);
vc_blank_period_ms = num2str(S_fig.blank_period_ms);
% P = get0_('P');
% nSites = numel(P.viSite2Chan);
csAns = inputdlg_({'blank_thresh (MAD)', 'blank_period_ms (millisecond)'}, ...
    'Common reference threshold', 1, {vc_blank_thresh, vc_blank_period_ms});
if isempty(csAns), return; end
% if isnan(site_start) || isnan(site_end), return; end;
blank_thresh = str2num(csAns{1});
blank_period_ms = str2num(csAns{2});
if isempty(blank_thresh), blank_thresh = 0; end
if isnan(blank_thresh) || isnan(blank_period_ms), return; end

S_fig = set_(S_fig, 'blank_thresh', blank_thresh, 'blank_period_ms', blank_period_ms);
set(hFig, 'UserData', S_fig);
S_fig = Fig_preview_update_(hFig, S_fig, 1);
end %func


%--------------------------------------------------------------------------
% 8/17/17 JJJ: created
function vrRef = mr2ref_(mnWav_filt, vcCommonRef, viSite_bad)
if nargin<2, vcCommonRef = 'mean'; end
if nargin<3, viSite_bad = []; end
if ~isempty(viSite_bad)
    mnWav_filt(:,viSite_bad) = [];
end
if strcmpi(vcCommonRef, 'median')
    vrRef = median(mnWav_filt, 2);
else
    vrRef = mean(mnWav_filt, 2);
end
end %func



%--------------------------------------------------------------------------
% 9/29/17 JJJ: delta==0 (found other spike that has exactly same PCA score) is set to nan
% 9/5/17 JJJ: Created
function [icl, x, z] = detrend_local_(S_clu, P, fLocal)
if nargin<3, fLocal = 1; end
maxCluPerSite = get_set_(P, 'maxCluPerSite', 20); % get 10 clu per site max
S0 = get0();
% cvi_rho_rank_site = cellfun(@(vi)rankorder_(S_clu.rho(vi), 'ascend'), S0.cviSpk_site, 'UniformOutput', 0);
% cvi_delta_rank_site = cellfun(@(vi)rankorder_(S_clu.delta(vi), 'ascend'), S0.cviSpk_site, 'UniformOutput', 0);

% detrend for each site and apply
if fLocal
    cvi_rho_site = cellfun(@(vi)S_clu.rho(vi), S0.cviSpk_site, 'UniformOutput', 0);
    cvi_delta_site = cellfun(@(vi)S_clu.delta(vi), S0.cviSpk_site, 'UniformOutput', 0);
%     [x0, y0] = deal(S_clu.rho, S_clu.delta);
    cvi_cl = cell(size(S0.cviSpk_site));
    x = log10(S_clu.rho);
    z = zeros(size(S_clu.delta), 'like', S_clu.delta);
    for iSite = 1:numel(S0.cviSpk_site)
        [rho_, delta_] = deal(cvi_rho_site{iSite}, cvi_delta_site{iSite});
        x_ = log10(rho_);
%         y_ = log10(rho_) + log10(delta_);
        y_ = (delta_);
        viDetrend = find(delta_ < 1 & delta_ > 0 & rho_ > 10^P.rho_cut & rho_ < .1 & isfinite(x_) & isfinite(y_));            
        [y_det_, z_] = detrend_(x_, y_, viDetrend, 1);
        viSpk_ = S0.cviSpk_site{iSite};        
        if isempty(viSpk_), continue; end
        vl_zero_ = find(delta_==0);
        [y_det_(vl_zero_), z_(vl_zero_)] = nan;
        [icl_, vrZ_] = find_topn_(y_det_, maxCluPerSite, ...
            find(rho_ > 10^P.rho_cut & ~isnan(y_det_)));
        if isempty(icl_), continue; end
        cvi_cl{iSite} = viSpk_(icl_); 
        z(viSpk_) = z_;
%         figure; plot(x_, y_det_, 'k.', x_(icl_), y_det_(icl_), 'ro'); axis([-5 0 -1 1]); grid on;
    end
    icl = cell2vec_(cvi_cl);
else
    x = log10(S_clu.rho);
    y = S_clu.delta;
    viDetrend = find(S_clu.delta < 1 & S_clu.delta > 0 & S_clu.rho > 10^P.rho_cut & S_clu.rho < .1 & isfinite(x) & isfinite(y));    
    [~, z] = detrend_(x, y, viDetrend, 1);
    z(S_clu.delta==0) = nan;
    [icl, vrZ_] = find_topn_(z, maxCluPerSite * numel(S0.cviSpk_site), ...
        find(S_clu.rho > 10^P.rho_cut & ~isnan(z)));
    icl(vrZ_ < 10^P.delta1_cut) = [];    
%     icl = find(x>=P.rho_cut & z>=10^P.delta1_cut);
end


if nargout==0
    figure; plot(x,z,'.', x(icl),z(icl),'ro'); grid on; 
    axis([-5 0 -20 100]);
    title(sprintf('%d clu', numel(icl)));
end
end %func


%--------------------------------------------------------------------------
function [y_, z] = detrend_(x, y, vi, fQuadratic)
% perform quadratic detrend and return z score
if nargin<4, fQuadratic = 1; end
if nargin<3, vi = 1:numel(x); end
x=x(:);
y=y(:);
% y = log10(y);
eps_ = eps(class(x));
if fQuadratic
    xx = [x.^2+eps_, x + eps_, ones(numel(x),1)];
else %linear detrending
    xx = [x + eps_, ones(numel(x),1)];
end
m = xx(vi,:) \ y(vi);
y_ = y - xx * m;
if nargout>=2
    z = (y_ - nanmean(y_(vi))) / nanstd(y_(vi));
end
end %func


%--------------------------------------------------------------------------
function  [viTop, vrTop] = find_topn_(vr, nMax, vi)
if nargin<3, vi = 1:numel(vr); end
nMax = min(nMax, numel(vi));
[~, viSrt] = sort(vr(vi), 'descend');
viTop = vi(viSrt(1:nMax));
vrTop= vr(viTop);
end %func


%--------------------------------------------------------------------------
function vr = cell2vec_(cvr)
% remove empty
cvr = cvr(:);
cvr = cvr(cellfun(@(x)~isempty(x), cvr));
for i=1:numel(cvr)
    vr_ = cvr{i};
    cvr{i} = vr_(:);
end
vr = cell2mat(cvr);
end %func


%--------------------------------------------------------------------------
% 09/14/17 JJJ: Export settings to a prm file 
function export_prm_(vcFile_prm, vcFile_out_prm, fShow)
% export current P structure to a parameter file
% export_prm_(vcFile_prm): read and write to the vcFile_prm (append default.prm)
% export_prm_(vcFile_prm, vcFile_out_prm): read from vcFile_prm and output to vcFile_out_prm
if nargin<3, fShow = 1; end

if isempty(vcFile_out_prm), vcFile_out_prm = vcFile_prm; end
copyfile(jrcpath_(read_cfg_('default_prm')), vcFile_out_prm, 'f');
P = get0_('P');
if isempty(P), P = file2struct_(vcFile_prm); end
edit_prm_file_(P, vcFile_out_prm);
vcMsg = sprintf('Full parameter settings are exported to %s', vcFile_out_prm);
fprintf('%s\n', vcMsg);
if fShow
    msgbox_(vcMsg);
    edit(vcFile_out_prm);
end
end %func


%--------------------------------------------------------------------------
% 9/17/17 JJJ: sample size is determined by the smallest file in the chan recording set
% 9/15/17 JJJ: Created
function vcFile_prm = import_intan_(vcFile_dat, vcFile_prb, vcArg3)

% vcFile_dat = 'E:\ZGao\Exp_20170710_1302_run000\amp-A-000.dat';
% vcFile_dat(end-7:end-5) = '000'; % replace to 000 index
csFiles_dat = sort(dir_file_(vcFile_dat, 1)); % sort by the channel name
[vcDir, vcFile] = fileparts(csFiles_dat{1});
vcFile_bin = [vcDir, '.bin'];
[~, vcFile_prb_] = fileparts(vcFile_prb);
vcFile_prm = strrep(vcFile_bin, '.bin', sprintf('_%s.prm', vcFile_prb_));
nChans = numel(csFiles_dat);
nBytes_file = min(cellfun(@(vc)getBytes_(vc), csFiles_dat));
P = struct('vcDataType', 'int16', 'probe_file', vcFile_prb, 'nChans', nChans, ...
    'uV_per_bit', .195, 'sRateHz', 30000, 'nBytes_file', nBytes_file);

nSamples = P.nBytes_file / bytesPerSample_(P.vcDataType); 

% Read file and output
mnWav  = zeros([nSamples, nChans], P.vcDataType);
for iFile = 1:numel(csFiles_dat)
    try
        fid_ = fopen(csFiles_dat{iFile}, 'r');
        mnWav(iFile:nChans:end) = fread(fid_, inf, ['*', P.vcDataType]);
        fclose(fid_);
        fprintf('Loaded %s\n', csFiles_dat{iFile});
    catch
        fprintf(2, 'error %s\n', csFiles_dat{iFile});
    end
end
write_bin_(vcFile_bin, mnWav);
clear mnWav;

% Write to a .prm file
try
    S_prb = file2struct_(vcFile_prb);
    if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
    if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
catch
    disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
end
P.duration_file = nSamples / P.sRateHz; %assuming int16
P.version = jrc_version_();
P.vcFile_prm = vcFile_prm;
P.vcFile = vcFile_bin;
copyfile(jrcpath_(read_cfg_('default_prm')), P.vcFile_prm, 'f');
edit_prm_file_(P, P.vcFile_prm);
vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
disp(vcPrompt);
edit(P.vcFile_prm);
end


%--------------------------------------------------------------------------
% 9/17/17 JJJ: Created for SPAARC
function vcFile_prm = import_nsx_(vcFile_nsx, vcFile_prb, vcTemplate_prm)
% Import neuroshare format
% sample size is determined by the smallest file in the chan recording set
if nargin<3, vcTemplate_prm = ''; end
if matchFileExt_(vcFile_prb, '.prm')
    vcTemplate_prm = vcFile_prb;
    S_ = file2struct_(vcTemplate_prm);
    vcFile_prb = S_.probe_file;
end

% vcFile_nsx = 'E:\TimBruns\Ichabod Trial 14\exp_9_ichabod0014.ns5';
if ~exist(vcFile_nsx, 'file'), error('File does not exist.') ;end
[P, nSamples] = nsx_info_(vcFile_nsx);
% [P, nSamples, vcFile_bin] = nsx2bin_(vcFile_nsx, 1);
% P.fInverse_file = 1;
% [mnWav, hFile, P] = load_nsx_(vcFile_nsx);
P.probe_file = vcFile_prb;
P.vcFile = vcFile_nsx;
% mnWav = mnWav * -1; %inverse polarity
[~, vcFile_prb_] = fileparts(vcFile_prb);
vcFile_prm = subsFileExt_(P.vcFile, sprintf('_%s.prm', vcFile_prb_));
if isempty(vcTemplate_prm)
    vcTemplate_prm = jrcpath_(read_cfg_('default_prm'));
end
assert(exist_file_(vcTemplate_prm), sprintf('Template file does not exist: %s', vcTemplate_prm));

% Write to a .prm file
try
    S_prb = file2struct_(find_prb_(vcFile_prb));
    if isfield(S_prb, 'maxSite'), P.maxSite = S_prb.maxSite; end
    if isfield(S_prb, 'nSites_ref'), P.nSites_ref = S_prb.nSites_ref; end
catch
    disperr_(sprintf('Error loading the probe file: %s\n', vcFile_prb));
end
P.duration_file = nSamples / P.sRateHz; %assuming int16
P.version = jrc_version_();
P.vcFile_prm = vcFile_prm;
% P.vcFile = vcFile_bin;
copyfile(vcTemplate_prm, P.vcFile_prm, 'f');
edit_prm_file_(P, P.vcFile_prm);
vcPrompt = sprintf('Created a new parameter file\n\t%s', P.vcFile_prm);
disp(vcPrompt);
edit(P.vcFile_prm);
end


%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function [mnWav, hFile, P] = load_nsx_(vcFile_nsx)
addpath('./neuroshare/');
[ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
% [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
nChans = sum(vlAnalog_chan);
% viElecID = double([hFile.Entity.ElectrodeID]);
P = struct('vcDataType', 'int16', 'nChans', nChans, ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sRateHz', 30000 / hFile.FileInfo.Period);

fprintf('Loading %s...', vcFile_nsx); t_load = tic;
fid = hFile.FileInfo.FileID;
if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
    fseek(fid, hFile.FileInfo.BytesHeaders + 9, -1);
else
    fseek(fid, hFile.FileInfo.BytesHeaders, -1);
end
mnWav = fread(fid, [nChans, nSamples], '*int16');  
fclose(fid);
fprintf('took %0.1fs\n', toc(t_load));

if nargout==0, assignWorkspace_(mnWav, hFile); end
end % func


%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function mrRate_clu = clu_rate_(S_clu, viClu, nSamples)
S0 = get0_();
P = S0.P;
if nargin<2, viClu = []; end 
if nargin<3, nSamples = []; end

if isempty(viClu), viClu = 1:S_clu.nClu; end
sRateHz_rate = get_set_(P, 'sRateHz_rate', 1000);
filter_sec_rate = get_set_(P, 'filter_sec_rate', 1);
nFilt = round(sRateHz_rate * filter_sec_rate / 2);   
filter_shape_rate = lower(get_set_(P, 'filter_shape_rate', 'triangle'));
switch filter_shape_rate
    case 'triangle'
        vrFilt = ([1:nFilt, nFilt-1:-1:1]'/nFilt*2/filter_sec_rate);
    case 'rectangle'
        vrFilt = (ones(nFilt*2, 1) / filter_sec_rate);
%     case 'interval3'
%         1/2 1/3 1/6. interval history
end %switch
vrFilt = single(vrFilt);

if isempty(nSamples), nSamples = round(P.duration_file * sRateHz_rate); end
mrRate_clu = zeros([nSamples, numel(viClu)], 'single');
for iClu1 = 1:numel(viClu)
    iClu = viClu(iClu1);
    viSpk_clu = S_clu.cviSpk_clu{iClu};
    viTime_clu = S0.viTime_spk(viSpk_clu);           
    viTime_ = round(double(viTime_clu) / P.sRateHz * sRateHz_rate);
    viTime_ = max(min(viTime_, nSamples), 1);
    mrRate_clu(viTime_, iClu1) = 1;        
    mrRate_clu(:,iClu1) = conv(mrRate_clu(:,iClu1), vrFilt, 'same');
end
end %func


%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function plot_aux_corr_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux, iCluPlot)
if nargin<5, iCluPlot = []; end
% show the firing rate and plot the 
[vrCorr_srt, viSrt] = sort(vrCorr_aux_clu, 'descend');
nClu = numel(vrCorr_aux_clu);
[P, S_clu] = get0_('P', 'S_clu');
P = loadParam_(P.vcFile_prm);
nClu_show = min(get_set_(P, 'nClu_show_aux', 4), nClu);
vcLabel_aux = get_set_(P, 'vcLabel_aux', 'aux');
nSubsample_aux = get_set_(P, 'nSubsample_aux', 100);
if ~isempty(iCluPlot)
   nClu_show = 1;
   viSrt = iCluPlot;
end
    
hFig = create_figure_('FigAux', [.5 0 .5 1], P.vcFile_prm,1,1);
hTabGroup = uitabgroup(hFig);
for iClu1 = 1:nClu_show
    iClu = viSrt(iClu1);
    htab1 = uitab(hTabGroup, 'Title', sprintf('Clu %d', iClu), 'BackgroundColor', 'w');    
    ax_ = axes('Parent', htab1);
    subplot(2, 1, 1);    
    ax_ = plotyy(vrTime_aux, mrRate_clu(:,iClu), vrTime_aux, vrWav_aux);
    xlabel('Time (s)');
    ylabel(ax_(1),'Firing Rate (Hz)');
    ylabel(ax_(2), vcLabel_aux);    
    iSite_ = S_clu.viSite_clu(iClu);
    vcTitle_ = sprintf('Clu %d (Site %d, Chan %d): Corr=%0.3f', ...
        iClu, iSite_, P.viSite2Chan(iSite_), vrCorr_aux_clu(iClu));
    title(vcTitle_);
    set(ax_, 'XLim', vrTime_aux([1,end]));
    grid on;
    
    subplot(2, 1, 2);
    plot(vrWav_aux(1:nSubsample_aux:end), mrRate_clu(1:nSubsample_aux:end,iClu), 'k.');
    xlabel(vcLabel_aux); 
    ylabel('Firing Rate (Hz)');
    grid on;
end %for
end %func


%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function plot_aux_rate_(fSelectedUnit)
% Aux channel vs. rate
if nargin<1, fSelectedUnit = 0; end %plot all
[P, S_clu, iCluCopy] = get0_('P', 'S_clu', 'iCluCopy');
P = loadParam_(P.vcFile_prm);
[vrWav_aux, vrTime_aux] = load_aux_(P);
if isempty(vrWav_aux), msgbox_('Aux input is not found'); return; end
mrRate_clu = clu_rate_(S_clu, [], numel(vrWav_aux));
vrCorr_aux_clu = arrayfun(@(i)corr(vrWav_aux, mrRate_clu(:,i), 'type', 'Pearson'), 1:size(mrRate_clu,2));
if ~fSelectedUnit, iCluCopy = []; end
plot_aux_corr_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux, iCluCopy);
vcMsg = assignWorkspace_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux);    
% msgbox_(vcMsg);
end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Imported from SPARC grant
function [vrWav_aux, vrTime_aux] = load_aux_(P)

[vrWav_aux, vrTime_aux] = deal([]);
if isempty(P.vcFile), msgbox_('Multi-file mode is currently not supported'); return; end
[~,~,vcExt] = fileparts(P.vcFile);
vcFile_aux = get_set_(P, 'vcFile_aux', '');
if isempty(vcFile_aux)
    switch lower(vcExt)
        case '.ns5', vcFile_aux = subsFileExt_(P.vcFile, '.ns2');
        case {'.bin', '.dat'}, vcFile_aux = P.vcFile;
        otherwise
            fprintf(2, 'Unable to determine the aux file. You must manually specify "vcFile_aux".\n');
            return;
    end
end
if ~exist_file_(vcFile_aux), return; end

[~,~,vcExt_aux] = fileparts(vcFile_aux);
switch lower(vcExt_aux)
    case '.ns2'
        iChan_aux = get_set_(P, 'iChan_aux', 1);
        [mnWav_aux, hFile_aux, S_aux] = load_nsx_(vcFile_aux);
        scale_aux = hFile_aux.Entity(iChan_aux).Scale * P.vrScale_aux;
        vrWav_aux = single(mnWav_aux(iChan_aux,:)') * scale_aux;        
        sRateHz_aux = S_aux.sRateHz;
    case '.mat'
        S_aux = load(vcFile_aux);
        csField_aux = fieldnames(S_aux);
        vrWav_aux = S_aux.(csField_aux{1});
        sRateHz_aux = get_set_(P, 'sRateHz_aux', P.sRateHz_rate);
    case {'.dat', '.bin'}
        iChan_aux = get_set_(P, 'iChan_aux', []);
        if isempty(iChan_aux), return; end
        mnWav_aux = load_bin_(vcFile_aux, P.vcDataType);
        vrWav_aux = single(mnWav_aux(iChan_aux:P.nChans:end)') * P.uV_per_bit * P.vrScale_aux;  
        sRateHz_aux = get_set_(P, 'sRateHz_aux', P.sRateHz);
    otherwise
        fprintf(2, 'vcFile_aux: unsupported file format: %s\n', vcExt_aux);
        return;
end %switch
if nargout>=2, vrTime_aux = single(1:numel(vrWav_aux))' / sRateHz_aux; end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Merged to Master. This function is not used, instead JRCLUST deals with the offset
% 9/22/17 JJJ: Created for SPARC. 
function [P, nSamples, vcFile_bin] = nsx2bin_(vcFile_nsx, fInvert)
if nargin<2, fInvert = 0; end

nBuffer = 1e8; % in bytes
addpath('./neuroshare/');
vcFile_bin = subsFileExt_(vcFile_nsx, '.bin');
[ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
% [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
nChans = sum(vlAnalog_chan);
% viElecID = double([hFile.Entity.ElectrodeID]);
P = struct('vcDataType', 'int16', 'nChans', nChans, ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sRateHz', 30000 / hFile.FileInfo.Period);
fprintf('Loading %s...', vcFile_nsx); t_load = tic;

fid = hFile.FileInfo.FileID;
fidw = fopen(vcFile_bin, 'w');

if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
    fseek(fid, hFile.FileInfo.BytesHeaders + 9, -1);
else
    fseek(fid, hFile.FileInfo.BytesHeaders, -1);
end
nLoad = ceil(nSamples*nChans / nBuffer);
for iLoad = 1:nLoad
    if iLoad==nLoad
        nBuffer_ = nSamples*nChans - (nLoad-1) * nBuffer;
    else
        nBuffer_ = nBuffer;
    end
    vnBuffer = fread(fid, nBuffer_, '*int16');
    if fInvert, vnBuffer = -vnBuffer; end
    fwrite(fidw, vnBuffer, 'int16');
end
fclose(fid);
fclose(fidw);
fprintf('took %0.1fs\n', toc(t_load));
end %func


%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function [fid, nBytes, header_offset] = fopen_nsx_(vcFile_nsx)
addpath('./neuroshare/');
[ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
% [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
nChans = sum(vlAnalog_chan);
nBytes = nSamples * nChans * 2;
fid = hFile.FileInfo.FileID;

if ismember(hFile.FileInfo.FileTypeID, {'NEURALCD', 'NEUCDFLT'})
    header_offset = hFile.FileInfo.BytesHeaders + 9;
else
    header_offset = hFile.FileInfo.BytesHeaders;
end
fseek(fid, header_offset, -1);
end %func


%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function [P, nSamples, hFile] = nsx_info_(vcFile_nsx)
addpath('./neuroshare/');
[ns_RESULT, hFile] = ns_OpenFile(vcFile_nsx);
vlAnalog_chan= strcmpi({hFile.Entity.EntityType}, 'Analog');
nSamples = hFile.TimeSpan / hFile.FileInfo.Period;
% viElecID = double([hFile.Entity.ElectrodeID]);
P = struct('vcDataType', 'int16', 'nChans', sum(vlAnalog_chan), ...
    'uV_per_bit', hFile.Entity(1).Scale, 'sRateHz', 30000 / hFile.FileInfo.Period);
end %func


%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function mrCor = chancor_(mn, P)
% mr = single(mn);
mrCor = zeros(size(mn), 'single');
nSites = numel(P.viSite2Chan);
nSites_spk = P.maxSite*2+1-P.nSites_ref;
miSites = P.miSites(2:nSites_spk,:);
for iSite = 1:nSites
    viSites_ = miSites(:,iSite);
    vrDist_ = pdist2(P.mrSiteXY(iSite,:), P.mrSiteXY(viSites_,:));
    vrDist_ = 1 ./ (vrDist_ / min(vrDist_));
    vrCor_ = bsxfun(@times, single(mn(:,iSite)), single(mn(:,viSites_))) * vrDist_(:);
    mrCor(:,iSite) = conv(vrCor_, [1,2,3,2,1]/3, 'same');
%     mrCor(:,iSrite) = vrCor_;
end
mrCor(mrCor<0)=0;
mrCor = int16(-sqrt(mrCor));
% figure; 
% subplot 121; imagesc(mrWav1(1:2000,:));
% subplot 122; imagesc(((mrCor1(1:2000,:))));
end %func


%--------------------------------------------------------------------------
% 9/29/17 JJJ: Displaying the version number of the program and what's used. #Tested
function [vcVer, vcDate, vcVer_used] = jrc_version_(vcFile_prm)
if nargin<1, vcFile_prm = ''; end
vcVer = 'v3.0.6';
vcDate = '10/8/2017';
vcVer_used = '';
if nargout==0
    fprintf('%s (%s) installed\n', vcVer, vcDate);
end
try 
    if isempty(vcFile_prm)
        P = get0_('P');
        if ~isempty(P)
            fprintf('\t%s used in %s\n', P.version, P.vcFile_prm);
            vcVer_used = P.version;
        end
    elseif exist_file_(vcFile_prm)
        P = loadParam_(vcFile_prm);
        fprintf('\t%s used in %s\n', P.version, vcFile_prm);
        vcVer_used = P.version;
    end
catch
    ;
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function git_pull_(vcVersion)
% https://github.com/drbenvincent/github-sync-matlab
% startDir = cd();
if nargin<1, vcVersion = ''; end

repoURL = 'https://github.com/jamesjun/JRCLUST';
% repoName = 'JRCLUST';
try
    if isempty(vcVersion)
        code = system('git pull');
    else
        code = system(sprintf('git reset --hard "%s"', vcVersion));
    end
catch
	code = -1;
end
if code ~= 0
    fprintf(2, 'Not a git repository. Please run the following command to clone from GitHub.\n');    
    fprintf(2, '\tsystem(''git clone %s.git myDest''\n', repoURL);
    fprintf(2, '\tReplace "myDest" with the desired installation location or omit to install in ./JRCLUST.\n', repoURL);
    fprintf(2, '\tYou may need to install git from https://git-scm.com/downloads.\n');  
else
    edit('change_log.txt');
end
end %func


%--------------------------------------------------------------------------
% 9/28/17 JJJ: Created and tested
function wiki_download_()
repoURL = 'https://github.com/jamesjun/JRCLUST.wiki.git';
repoName = 'wiki';
if isempty(dir(repoName))
    vcEval = sprintf('git clone %s %s', repoURL, repoName);
else
    vcEval = sprintf('git pull %s', repoName);
end
try
	code = system(vcEval);
catch
	code = -1;
end
if code == 0
    fprintf('Wiki on GitHub is downloaded to ./%s/\n', repoName);
else
    fprintf(2, 'Please install git from https://git-scm.com/downloads.\n');        
end
end %func


%--------------------------------------------------------------------------
% 9/26/17 JJJ: Created and tested
function flag = exist_file_(vcFile)
% Different from exist(vcFile, 'file') which uses search path
if isempty(vcFile)
    flag = 0; 
else
    flag = ~isempty(dir(vcFile));
end
end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function gui_(vcArg1, vcFile_prm_)
% JRCLUST GUI interface

if ~isempty(vcArg1)
    vcFile_prm = vcArg1;
elseif ~isempty(vcFile_prm_)
    vcFile_prm = vcFile_prm_;
else
    vcFile_prm = '';
end
S_gui.vcFile_prm = vcFile_prm;
set0_(S_gui);
jrc3_gui();

end %func


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function wiki_(vcPage)
if nargin<1, vcPage = ''; end
if isempty(vcPage)
    web('https://github.com/jamesjun/JRCLUST/wiki'); 
else
    web(['https://github.com/jamesjun/JRCLUST/wiki/', vcPage]); 
end
end


%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created
function issue_(vcMode)
switch lower(vcMode)
    case 'post', web('https://github.com/jamesjun/JRCLUST/issues/new')
    case 'search', web('https://github.com/jamesjun/JRCLUST/issues')
end %switch
end %func


%--------------------------------------------------------------------------
% 10/8/17 JJJ: Created
function edit_(vcFile)
% vcFile0 = vcFile;
if ~exist_file_(vcFile)
    if matchFileExt_(vcFile, '.prb')
        vcFile = find_prb_(vcFile);
    else
        vcFile = jrcpath_(vcFile, 1);
    end
end
fprintf('Editing %s\n', vcFile);
edit(vcFile);
end %func
