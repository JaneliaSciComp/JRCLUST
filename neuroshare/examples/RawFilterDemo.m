% MATLAB Native Neuroshare API demo. Opens a neural recording and 
% demonstrates the application of a low pass filter (using Matlab's
% butter and filter functions) and downsampling the data (using
% Matlab's resample function).  This demo will take a file that contains
% raw and LFP data and by applying the filter and down sampling the raw
% data will reproduce the LFP.
%
% This demo should run well as along as it is run on a recording that 
% contains both LFP and raw continuous recordings.
% 
% For more information on the Neuroshare API see:
% http://neuroshare.sourceforge.net/index.shtml
%
% Ripple LLC 2013
% Elliott Barcikowski

%
% This functin uses MATLAB cell mode. Matlab cells are paritioned with the 
% double percent sign (%%). You can execute the code in a cell by moving 
% the cursor to the cell and pressing CTRL-Enter
%
% Contact: support@rppl.com
%% Neuroshare Exists
% Check to make sure neuroshare functions are in the MATLAB path
if ~exist('ns_CloseFile.m', 'file')
    disp('WARNING => the neuroshare functions were not found in your MATLAB path.');
    disp('looking in directory ../');    
    path(path, '../');    
    
    if exist('ns_CloseFile.m', 'file')
        disp('This script found the neuroshare functions in ../');
    else        
        disp('In a standard Trellis installation these functions are located in:');
        disp('   <trellis_install_dir>/Tools/matlab/neuroshare>');
    end
    
    % add the path as it is in SVN
    path(path, '..');
    
    if ~exist('ns_CloseFile.m', 'file')
        return;
    end
end

%% ns_CloseFile
% when finished using neuroshare functions or to open another data set
% first run ns_CloseFile(hFile). also it is good practice to clear the
% hFile variable from the workspace.
if exist('hFile', 'var')
    ns_CloseFile(hFile);
end

%% get file handle (hFile)
clear all
close all

% Get screen dimensions so that the sizes of the figure windows can
% be controlled for any screen resolution
SCREEN_SIZE = get(0, 'ScreenSize');
% modifiying header and footer parameters may beautify window sizes
SCREEN_FOOTER = 55;
SCREEN_HEADER = 150;
SCREEN_L_MARGIN = 5;
SCREEN_R_MARGIN = 5;

%% Open a neural recording
[rc, hFile] = ns_OpenFile;

[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);

% check to make sure the operation succeeded
if strcmp(ns_RESULT, 'ns_OK') ~= 1
    disp(['ERROR: ns_OpenFile() returned ' ns_RESULT]);
    return;
end

%% Get entity info
nsEntityInfo(nsFileInfo.EntityCount,1).EntityLabel = '';
nsEntityInfo(nsFileInfo.EntityCount,1).EntityType  = 0;
nsEntityInfo(nsFileInfo.EntityCount,1).ItemCount   = 0;
% The entity inforamtion is read using ns_GetEntityInfo()
for i = 1:nsFileInfo.EntityCount
    [~, nsEntityInfo(i,1)] = ns_GetEntityInfo(hFile, i);
end

%% Find entities with continuous data.
AnalogEntityIDs = find([nsEntityInfo.EntityType]==2);
ElecIDsWithContData =  unique([hFile.Entity(AnalogEntityIDs).ElectrodeID]);

% pick a random electrode that has continuous data
WantedElec = ElecIDsWithContData(randi(length(ElecIDsWithContData)));
WantedEntityIDs = find([hFile.Entity.ElectrodeID]==WantedElec);

lfpEntityID = WantedEntityIDs(2);
rawEntityID = WantedEntityIDs(3);

lfpCount = nsEntityInfo(lfpEntityID).ItemCount;
rawCount = nsEntityInfo(rawEntityID).ItemCount;

lfpTimes = (1:lfpCount);
rawTimes = (1:rawCount)/30;

[~, ~, lfpData] = ns_GetAnalogData(hFile, lfpEntityID, 1, lfpCount);
[~, ~, rawData] = ns_GetAnalogData(hFile, rawEntityID, 1, rawCount);

%% Apply low pass filter and downsampling to raw data to reproduce LFP.
% The filter parameters shown here will reproduce the filtering that
% occurs before the LFP data is recorded.

Nq = 15000; % Nyquist frequency => sampling rate / 2 

% Get filter coeffients
[zLow, pLow] = butter(4, 250/Nq, 'low');
% Apply filter
filtData = filter(zLow, pLow, rawData);
% Downsample the data.
simLFP = resample(filtData, 1, 30);


%% Draw Orginal Raw and LFP waveforms
figOriginal = figure('Name', 'Recorded Raw and LFP',...
    'Position', [ SCREEN_L_MARGIN SCREEN_FOOTER+SCREEN_SIZE(4)/4 ...
    SCREEN_SIZE(3)/3 SCREEN_SIZE(4)*2/3 - SCREEN_HEADER ], ...
    'NumberTitle', 'off');
    

subplot(2, 1, 1);      
plot(rawTimes, rawData);
title('Recorded Raw Data - \DeltaT = 30 kHz');
xlabel('[ms]');
ylabel('[\muV]');
subplot(2, 1, 2);
plot(lfpTimes, lfpData);
title('Recorded LFP Data - Low Pass at 250 Hz, \DeltaT = 1kHz');
xlabel('[ms]');
ylabel('[\muV]');


%% Draw the filtered and downsampled waveforms
figFiltered = figure('Name', 'Recorded Raw and LFP',...
    'Position', [ SCREEN_L_MARGIN+SCREEN_SIZE(3)/3 ...
    SCREEN_FOOTER+SCREEN_SIZE(4)/4 SCREEN_SIZE(3)/3 ...
    SCREEN_SIZE(4)*2/3 - SCREEN_HEADER], ...
    'NumberTitle', 'off');

xrange = [min(lfpTimes) max(lfpTimes)];

subplot(3, 1, 1);
plot(lfpTimes, simLFP);

set(gca, 'XLim', xrange);
title('Filtered and Downsampled Raw Data - Low Pass at 250 Hz and \DeltaT = 1kHz');
xlabel('[ms]');
ylabel('[\muV]');

subplot(3, 1, 2);
plot(lfpTimes, lfpData);
set(gca, 'XLim', xrange);
title('Recorded LFP Data - Low Pass at 250 Hz, \DeltaT = 1kHz');
xlabel('[ms]');
ylabel('[\muV]');

subplot(3, 1, 3);
plot(lfpTimes, simLFP-lfpData);
set(gca, 'XLim', xrange);
set(gca, 'YLim', [-1, 1]);
title('Difference between LFP and Filtered Raw Data');
xlabel('[ms]');
ylabel('[\muV]');


