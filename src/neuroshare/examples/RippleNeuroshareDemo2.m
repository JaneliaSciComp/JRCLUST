% Neuroshare matlab demo script
%
% Ripple LLC 2012
% support@rppl.com
% this demo uses cell mode. Matlab Editor Cells are paritioned using a
% double percent sign at the beging of the line (%%). to run the code in a
% particular cell you simply press ctrl-Enter
% processing variables that are useful for all sections
% The following are useful for orienting figure windows

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

[ns_RESULT, hFile] = ns_OpenFile;
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, 'Failed to open NEV or NSX data\n');
    return
end

%%  
% get file info structure. This is necessary particularly for getting the
% EntityCount or the number of Entities in the data set. After getting the
% number of Entities, we get all of the EntityInfo Structures for each
% Entity. The EntityInfo Structures contains the EntityType, EntityLabel,
% and the ItemCount. The ItemCount is the number of occurances or samples
% of each Entity. 
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
% get entity info structure. In order to access the EntityInfo
% information quickly is is written into a structure. To create the
% structure it must be preallocated first:
nsEntityInfo(nsFileInfo.EntityCount,1).EntityLabel = '';
nsEntityInfo(nsFileInfo.EntityCount,1).EntityType = 0;
nsEntityInfo(nsFileInfo.EntityCount,1).ItemCount = 0;
% the structure is then filled using ns_GetEntityInfo
for i = 1:nsFileInfo.EntityCount
    [~, nsEntityInfo(i,1)] = ns_GetEntityInfo(hFile, i);
end

%% get Event/Analog/Segment EntityIDs
% This gives the EntityID numbers for the Event, Analog, and Segment data.
% According to the neuroshare specificaion The EntityType is defined as
% follows:
%   
%                                  Unknown entity                   0
%   
%                                  Event entity                     1
%   
%                                  Analog entity                    2
%   
%                                  Segment entity                   3
%   
%                                  Neural event entity              4
% In the Ripple Specification for the neuroshare functions for nev/nsx
% files there are no Unknown Entities.  In the future. putting brackets 
% around a structure element concatinates it to an array: [struct.element1] 
% is an array of all of the element1 of the structure. Setting the array 
% to equal to a number creates a logical index of all of the elements 
% that equal the particular number. using the find function on the logical 
% index array gives a numerical index which intern is the EntityID of the 
% particular EntityType.
EventEntityID = find([nsEntityInfo.EntityType]==1);
AnalogEntityID = find([nsEntityInfo.EntityType]==2);
SegmentEntityID = find([nsEntityInfo.EntityType]==3);

% IJM: EntityLabel has changed to user defined/map file defined values.
% Use FileType to get 30 kS/s streams and ensure ElectrodeID does not
% include AIO
%find(~cellfun('isempty', strfind({nsEntityInfo.EntityLabel}, '30 kS/s')));
RawEntityID = find([hFile.Entity(:).FileType] == 4);
RawEntityID = ...
    RawEntityID([hFile.Entity(RawEntityID).ElectrodeID] < 513);

% IJM: EntityLabel has changed to user defined/map file defined values.
% Use FileType to get LP streams and ensure ElectrodeID does not
% include AIO
%find(~cellfun('isempty', strfind({nsEntityInfo.EntityLabel}, '30 kS/s')));
LFPEntityID = find([hFile.Entity(:).FileType] == 3);
LFPEntityID = ...
   LFPEntityID([hFile.Entity(LFPEntityID).ElectrodeID] < 513);

%% get Event/Analog/Segment EntityLabels
% This organizes the EntityLabels found in nsEntityInfo by the Entity type 
% for plot labeling purposes
EventLabel   = {nsEntityInfo(EventEntityID).EntityLabel};
AnalogLabel  = {nsEntityInfo(AnalogEntityID).EntityLabel};
SegmentLabel = {nsEntityInfo(SegmentEntityID).EntityLabel};

%% draw spike waveforms
% setup figure to draw spike waveforms for the first few electrodes (up
% to the first four).  Colors the waveforms based on spike sorting
% classifications.

% Pack the wanted positions into an array.  This figure will take up
% the left 1/3 of the width of the screen
position = [ SCREEN_L_MARGIN SCREEN_FOOTER ...
    SCREEN_SIZE(3)/3 SCREEN_SIZE(4) - SCREEN_HEADER ];
% if this script has been run before close the figure 
% associated with this section.
if exist('spikeHandle', 'var')
    close(spikeHandle)
end
% open spike figure
spikeHandle = figure('Position', position, 'NumberTitle', 'off', ...
    'Name', 'Spike Waveforms');

% use up to the first four electrodes found in the nev data
nplots = min(4, length(SegmentEntityID));
for iSegment=1:nplots
    subplot(nplots, 1, iSegment);
    % the details of plotting a given electrode spikes are packed in
    % this subfunction
    plot_spikes(hFile, SegmentEntityID(iSegment));
end

% produce a title textbox centered at the top of the center frame 
annotation(spikeHandle,'textbox',...
    [0.38 0.960 0.35 0.035],...
    'String',{'Spike Waveforms'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

%% draw spike time histogram
% Draw the spike times as a histogram for the first few electrodes (up
% to the first four).

% Pack the wanted positions into an array.  This figure will take up
% the center 1/3 of the width of the screen.
position = [ SCREEN_L_MARGIN+SCREEN_SIZE(3)/3 SCREEN_FOOTER ...
    SCREEN_SIZE(3)/3 SCREEN_SIZE(4) - SCREEN_HEADER ];
% if this script has been run before close the figure 
% associated with this section.
if exist('histHandle', 'var')
    close(histHandle)
end
% open spike histogram figure
histHandle = figure('Position', position, 'NumberTitle', 'off', ...
    'Name', 'Spike Times');

% use up to the first four electrodes found in the nev data
nplots = min(4, length(SegmentEntityID));
for iSegment=1:nplots
    subplot(nplots, 1, iSegment);
    % the details of plotting and producing the histograms for
    % given electrode spikes are packed in this subfunction
    spike_hist(hFile, SegmentEntityID(iSegment), 1);
end

% produce a title textbox centered at the top of the center frame 
annotation(histHandle,'textbox',...
    [0.40 0.960 0.35 0.035],...
    'String',{'Spike Times'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');



%% draw ISI histogram
% [ELB-2013/03/29] Three sets of plots is enough.  Removing this for the
% time being.  It can always be put back in later.
% % Draw histograms of the ISI(Inter-Spike Interval)times for the first
% % few electrodes (up to four).
% 
% % Pack the wanted positions into an array.  This figure will take up
% % the center 1/3 of the width of the screen.  This figure will end
% % up on top of the spike times histogram.
% position = [ SCREEN_L_MARGIN+SCREEN_SIZE(3)/3 SCREEN_FOOTER ...
%     SCREEN_SIZE(3)/3 SCREEN_SIZE(4) - SCREEN_HEADER ];
% % if this script has been run before close the figure 
% % associated with this section.
% if exist('isiHandle', 'var')
%     close(histHandle)
% end
% % open the isi figure  
% isiHandle = figure('Position', position, 'NumberTitle', 'off', ...
%     'Name', 'Spike Times');
% % use up to the first four electrodes found in the nev data
% nplots = min(4, length(SegmentEntityID));
% for iSegment=1:nplots
%     subplot(nplots, 1, iSegment);
%     % the details of plotting and producing the ISI histograms are
%     % found in the plot_isi_hist sub-function
%     plot_isi_hist(hFile, SegmentEntityID(iSegment));
% end
% 
% % produce a title textbox centered at the top of the center frame 
% annotation(isiHandle,'textbox',...
%     [0.4 0.960 0.35 0.035],...
%     'String',{'ISI Histogram'},...
%     'FontSize',14,...
%     'FitBoxToText','off',...
%     'LineStyle','none');

%% draw analog data
% Plot the continuous data for the first two electrodes.  We plot both the
% LFP and raw (if they are present).

plotPeriodSec = 20; % Amount of data to plot in seconds
rawPoints = 20*30000; % number of raw data points that we need at 30kHz
lfpPoints = 20*1000; % number of LFP points at 1kHz
% Pack the wanted positions into an array.  This figure will take up
% the right 1/3 of the width of the screen.
    position = [ SCREEN_L_MARGIN+2*SCREEN_SIZE(3)/3 SCREEN_FOOTER ...
    SCREEN_SIZE(3)/3 SCREEN_SIZE(4) - SCREEN_HEADER ];
% if this script has been run before close the figure 
% associated with this section.
if exist('analogHandle', 'var')
    close(analogHandle)
end
% open the analog figure  
analogHandle = figure('Position', position, 'NumberTitle', 'off', ...
    'Name', 'Continuous Waveforms');

% In the first pannel plot the full analog 
% waveform for the first analog channel
subplot(4, 1, 1)
plot_analog(hFile, RawEntityID(1), 1,...
    min(nsEntityInfo(RawEntityID(1)).ItemCount, rawPoints));
% In the second panel plot the first 500 
% points of the first analog channel
subplot(4, 1, 2)
plot_analog(hFile, LFPEntityID(1), 1, ...
    min(nsEntityInfo(LFPEntityID(1)).ItemCount, lfpPoints));

subplot(4, 1, 3)
plot_analog(hFile, RawEntityID(2), 1,...
    min(nsEntityInfo(RawEntityID(2)).ItemCount, rawPoints));

subplot(4, 1, 4)
plot_analog(hFile, LFPEntityID(2), 1, ...
        min(nsEntityInfo(LFPEntityID(2)).ItemCount, lfpPoints));

% produce a title textbox centered at the top of the center frame 
annotation(analogHandle,'textbox',...
    [0.30 0.960 0.55 0.035],...
    'String',{'Continuous Waveforms'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');
