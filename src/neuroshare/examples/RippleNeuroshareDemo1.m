% MATLAB Native Neuroshare API demo. Opens a neural recording and 
% demonstrates various examples of how to access and visualize data
%
% For more information on the Neuroshare API see:
% http://neuroshare.sourceforge.net/index.shtml
%
% Ripple LLC 2011-2012
% Mohammad Raza
% Andrew Wilder
%
% This functin uses MATLAB cell mode. Matlab cells are paritioned with the 
% double percent sign (%%). You can execute the code in a cell by moving 
% the cursor to the cell and pressing CTRL-Enter
%
% Contact: support@rppl.com

function RippleNeuroshareDemo1(varargin)
% neuroshareDemo - Opens and plots data from a recording
%
% Usage:
% neuroshareDemo('recordingPath')
%
% Parameters:
% recordingPath  A string that specifies the path to a recording. The 
%                recording can be a specific file or the base name of a set
%                of files that all belong to the same recording.
%
% Remarks: 
% Neural data from a single recording session may be stored in multiple
% files. In NEV recordings, for example, events and spikes are stored in a 
% file with a .nev extension while continuous, analog data are stored in 
% files with .ns* extensions where the * corresponds to the sampling 
% period. The Neuroshare API is agnostic to file structure. In the MATLAB 
% Native Neuroshare implementation, the default assmption is that a  
% recording is split across multiple files. Thus when an individual file
% is passed as an argument, the ns_OpenFile() method checks the directory  
% for other files from the same recording. If such files are found they are 
% also opened. (NOTE: Files with the same base name are assumed to be from 
% the same recording).
%

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


%% Create a neuroshare file handle

% This is a necessary step; the file handle is used as an argument in all 
% subsequent neuroshre functions.
if nargin > 0
    [ns_RESULT, hFile] = ns_OpenFile(varargin{1});
else
    [ns_RESULT, hFile] = ns_OpenFile();
end

% check to make sure the operation succeeded
if( strcmp(ns_RESULT, 'ns_OK') ~= 1 )
    disp(['ERROR: ns_OpenFile() returned ' ns_RESULT]);
    return;
end

%% Read high-level information about the recording

% The file info structure "nsFileInfo" contains information about when the 
% recording was created, the timespan of the recording, timestamp 
% resolution in seconds, and the number of entities in the recording. 
[ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);

% check to make sure the operation succeeded
if( strcmp(ns_RESULT, 'ns_OK') ~= 1 )
    disp(['ERROR: ns_GetFileInfo() returned ' ns_RESULT]);
    return;
end

% Using the entity count from the nsFileInfo structure we can preallocate 
% memory for storing the nsEntityInfo structure for each entity. The 
% nsEntityInfo structure contains the entity label, type, and item count 
% (NOTE: Item count is the number of events or samples recorded for an 
% entity).
nsEntityInfo(nsFileInfo.EntityCount,1).EntityLabel = '';
nsEntityInfo(nsFileInfo.EntityCount,1).EntityType  = 0;
nsEntityInfo(nsFileInfo.EntityCount,1).ItemCount   = 0;

% The entity inforamtion is read using ns_GetEntityInfo()
for i = 1:nsFileInfo.EntityCount
    [~, nsEntityInfo(i,1)] = ns_GetEntityInfo(hFile, i);
end


%% Organize all entities in the recording by type

% As per the Neuroshare specification, The codes for the various types of 
% entities are as follows:
%   
%               Unknown entity                   0
%   
%               Event entity                     1
%   
%               Analog entity                    2
%   
%               Segment entity                   3
%   
%               Neural Event entity              4
%
% NOTES: 
% 1) In the data files saved by Trellis there currently are no entities
% of type Unknown or Neural Event. Neural Event type entities will be 
% present when online sorting is enabled (in a future Grapevine release)
%
% 2) A review of MATLAB syntax:
%
% - Enclosing a structure element in brackets concatinates it to form an
%   array. For example, [struct.element1] is an array of all of the  
%   element1 of the structure struct. 
%
% - Setting an array equal to a number creates a logical index of all of  
%   the elements that equal the particular number. 
% 
% - Using the find function on the logical index array gives a numerical 
%   index which is the EntityID of the particular EntityType.
EventEntityIDs   = find([nsEntityInfo.EntityType]==1);
AnalogEntityIDs  = find([nsEntityInfo.EntityType]==2);
SegmentEntityIDs = find([nsEntityInfo.EntityType]==3);


%% Organize entity labels by entity type

% Labels are usefull for identifying entities and labeling plots.
EventEntityLabels   = {nsEntityInfo(EventEntityIDs).EntityLabel};
AnalogEntityLabels  = {nsEntityInfo(AnalogEntityIDs).EntityLabel};
SegmentEntityLabels = {nsEntityInfo(SegmentEntityIDs).EntityLabel};


%% Organize item counts by entity type

% This is usefull for any data retrieval using the ns_Get*Data functions. 
% For example, one could access data for the third Analog Entity with the
% following line:
%
% [ns_RESULT, SampleCount, Samples] = ...
%      ns_GetAnalogData(hFile, AnalogEntityIDs(3), 1, AnalogItemCounts(3));
%
EventItemCounts   = [nsEntityInfo(EventEntityIDs).ItemCount];
AnalogItemCounts  = [nsEntityInfo(AnalogEntityIDs).ItemCount];
SegmentItemCounts = [nsEntityInfo(SegmentEntityIDs).ItemCount];


%% Get data and timestamps for Event-type Entities

% Before reading the event data and their corresponding time stamps we 
% first pre-allocate arrays to hold all of the information (time stamps, 
% data sizes, and the data itself). We then loop throuh the list of event-
% type entities populating these arrays.
EventEntityCount  = length(EventEntityIDs);
Events            = nan(max(EventItemCounts), EventEntityCount);
EventTimes        = nan(max(EventItemCounts), EventEntityCount);
EventSizes        = nan(max(EventItemCounts), EventEntityCount);
for i = 1:EventEntityCount
    for j = 1:EventItemCounts(i)
        [~, EventTimes(j,i), Events(j,i), EventSizes(j,i)] = ...
            ns_GetEventData(hFile, EventEntityIDs(i), j);
    end
end


%% Get waveforms and sample counts for Analog-type entities

% When the number of Analog samples recorded is large it is possible to 
% exceed the memory limit of your computer if you try to load all samples 
% at once. Here we demonstrate how to limit the number of samples accessed 
% for a specific entity. Use MaxAnalogSampleCount to control the size of
% the sample set.
AnalogEntityCount = length(AnalogEntityIDs);
MaxAnalogSampleCount = 10*30000; % 10 seconds @ 30,000 samp/sec
AnalogSampleReadStart = 1;
% Calculate the number of samples to read for each analog entity.  Will
% hold the max number of data to be read (MaxAnalogSampleCount) or the
% number of samples in the data, which ever is smaller.
AnalogSampleReadCounts = min(AnalogItemCounts, MaxAnalogSampleCount);
% Analog entities may have various item counts (sample length) depending on
% the rate of sampling.  As the data is retrieved, arrays are allocated for
% the fastest sampling (largest number of samples).  The data correspoding 
% to slower sampling will be padded with zeros.
AnalogSampleCounts = zeros(AnalogEntityCount, 1);
AnalogWaveforms = zeros(max(AnalogSampleReadCounts), AnalogEntityCount);
fprintf('retrieving continuous data for channel:\n');
for i = 1:AnalogEntityCount
    fprintf('%3d, ', i);
    if mod(i, 20) == 0
        fprintf('\n');
    end
    % Read data from files and pack it in the 
    [~, AnalogSampleCounts(i), AnalogWaveforms(1:AnalogSampleReadCounts(i),i)] = ...
        ns_GetAnalogData(hFile, AnalogEntityIDs(i), ...
                         AnalogSampleReadStart, AnalogSampleReadCounts(i));
end
fprintf('\n');
% NOTE: The NEV specification supports pauses in a recording. The 
% Neuroshare API does not make specific provisions for handline pauses. 
% Thus for NEV continuous reordings with pauses, the MATLAB Native 
% Neuroshare implementation of ns_GetAnalogData() returns the number of 
% samples before the first pauses.


%% Get waveforms and timestamps for all Segment-type entities

SegmentEntityCount   = length(SegmentEntityIDs);
SegmentTimeStamps    = nan(SegmentEntityCount, max(SegmentItemCounts));
SegmentUnitIDs       = nan(SegmentEntityCount, max(SegmentItemCounts));
SegmentSampleCounts  = nan(SegmentEntityCount, max(SegmentItemCounts));
SegmentWaveforms     = cell(size(SegmentTimeStamps));
fprintf('retrieving spikes for channel:\n')
for i = 1:SegmentEntityCount
    fprintf('%3d, ', i);
    if mod(i, 20) == 0
        fprintf('\n');
    end
    for j = 1:SegmentItemCounts(i)
        [~, ...
         SegmentTimeStamps(i,j), ...
         SegmentWaveforms{i,j}, ...
         SegmentSampleCounts(i,j), ...
         SegmentUnitIDs(i,j)] = ...
             ns_GetSegmentData(hFile, SegmentEntityIDs(i), j);
    end
end
fprintf('\n');

%% Examples of data visualization
scrnsize = get( 0, 'ScreenSize' );

% set up the sizes of the figures
figAnalogPos      = [(scrnsize(3)*0.02)/2, scrnsize(4)*0.8-90, ...
                      scrnsize(3)*0.98, scrnsize(4)*0.2];

figSegmentPos     = [figAnalogPos(1), 60, ...
                      scrnsize(3)*0.35, figAnalogPos(2)-160];
             
figSegmentGridPos = ...
            [figSegmentPos(1) + figSegmentPos(3) + 20, ...
             figSegmentPos(2), ...
             (figAnalogPos(1)+figAnalogPos(3))-(figSegmentPos(1)+ ...
             figSegmentPos(3)+20), ...
             figSegmentPos(4)];
             
% -           
% - 1) Plot some Analog waveforms
% -

% make sure at least one analog entity exists
if AnalogEntityCount > 0
    
    % create the figure
    figAnalog = figure('Name', ...
                       'Example plot of Analog Entity Recording', ...
                       'NumberTitle', ...
                       'off');
    set(figAnalog, 'Position', figAnalogPos);
    analogEntityIdx       = 1;
    analogEntityID        = AnalogEntityIDs(analogEntityIdx);
    analogEntityLabel     = AnalogEntityLabels{analogEntityIdx};
    analogEntityWaveform  = AnalogWaveforms(:,analogEntityIdx);
    analogDataUnit        = 'uV';
    analogTimeUnit        = 'ms';
    sampPerAnalogTimeUnit = 30;
    analogTimeVect        = ...
        ((0:(length(analogEntityWaveform)-1))+(AnalogSampleReadStart-1))./sampPerAnalogTimeUnit;

    % plot the waveform
    plot(analogTimeVect*1000, analogEntityWaveform);
    xlabel(['time (' analogTimeUnit ')']);
    ylabel(analogDataUnit);
    title(['Data for Entity [' analogEntityLabel ']']);
    xlim([analogTimeVect(1), analogTimeVect(end)]);
    
end

% -
% - 2) Plot segment waveforms and activity histogram for a specific entity.
% -

% make sure at least one Segment  entity exists
if SegmentEntityCount > 0

    figSegment = figure('Name', ...
                        'Example plot of Segment Entity Recording', ...
                        'NumberTitle', ...
                        'off');
    set(figSegment, 'Position', figSegmentPos);

    % In NEV recordings, all of the spike waveforms are of the same length.
    % Here we show an example of how to plot all of the spike waveforms for a 
    % specific segment-type entity.
    entityIdx = 1;
    entityID = SegmentEntityIDs(entityIdx);
    entityLabel = SegmentEntityLabels{entityIdx};
    spikes = SegmentWaveforms(entityIdx,:);
    spikeCount = SegmentItemCounts(entityIdx);
    spikeTimes = SegmentTimeStamps(entityIdx,:);
    spikeTimes = spikeTimes(~isnan(spikeTimes)); % remove NaN values
    sampCount = SegmentSampleCounts(entityIdx,1);
    spikeUnit = 'uV';
    spikeTimeUnit = 'ms';
    sampPerTimeUnit = 30; % NEV spike sample rate is 30,000 samp/sec
    timeBaseVect = (0:(sampCount-1))./sampPerTimeUnit;
    histTimeUnit = 'sec';
    histBinCenters = (((1:10)-.5)./10)*spikeTimes(end);

    % cell array that holds the time vector for each spike waveform.
    timeVects = cell(size(spikes));
    [timeVects{1,1:spikeCount}] = deal(timeBaseVect);

    % plot all the waveforms
    subplot(2, 1, 1);
    hold on
    cellfun(@plot, timeVects, spikes)
    title(['Spike events for Entity [' entityLabel ']']);
    xlabel(['Time (' spikeTimeUnit ')']);
    ylabel(spikeUnit);
    xlim([timeBaseVect(1), timeBaseVect(end)]);

    % plot the histogram
    subplot(2, 1, 2);
    hist(spikeTimes, histBinCenters)
    title(['Spike time histogram for Entity [' entityLabel ']']);
    xlabel(['Time (' histTimeUnit ')']);
    ylabel('count');
    xlim([0, histBinCenters(end)+histBinCenters(1)]);
    
end
    
    
% -
% - 3) Spike waveforms for all entities
% -

% make sure at least one Segment  entity exists
if SegmentEntityCount > 0
    
    % make the figure
    figSegmentGrid = ...
      figure('Name', ...
             ['Example plot of several' ...
             ' Spikes for all Segment Entities'], ...
             'NumberTitle', ...
             'off');
    set(figSegmentGrid, 'Position', figSegmentGridPos);
    
    rowCount = ceil(sqrt(SegmentEntityCount));
    colCount = ceil(SegmentEntityCount/rowCount);
    
    % plot the first spikeCount spikes for each entity
    for i=1:SegmentEntityCount
        
        entityIdx        = i;
        entityLabel      = SegmentEntityLabels{entityIdx};
        spikes           = SegmentWaveforms(entityIdx,:);
        spikesValidFlg   = ~cellfun(@isempty, spikes);
        lastValidIdx     = find(spikesValidFlg==0, 1, 'first')-1;
        
        % catch case where all spikes are valid
        if isempty(lastValidIdx)
            lastValidIdx = length(spikes);
        end
        
        % change this to plot more/fewer spikes
        spikeCount       = min(10, lastValidIdx); 
        spikes           = spikes(1:spikeCount);
        sampCount        = SegmentSampleCounts(entityIdx,1);
        sampPerMSec      = 30; % NEV spike sample rate is 30,000 samp/sec
        timeBaseVect     = (0:(sampCount-1))./sampPerMSec;

        % cell array that holds the time vector for each waveform.
        timeVects = cell(size(spikes));
        [timeVects{1,1:spikeCount}] = deal(timeBaseVect);

        % plot the selected waveforms
        subplot(rowCount, colCount, i);
        hold on
        cellfun(@plot, timeVects, spikes)
        title(entityLabel);
        xlim([timeBaseVect(1), timeBaseVect(end)]);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        
    end

end

%% Close the file handle

% When you are finished using neuroshare functions or before you open 
% another recording you must call ns_CloseFile(). It is also good practice 
% to clear the hFile variable from the workspace.
ns_RESULT = ns_CloseFile(hFile);

% check to make sure the operation succeeded
if( strcmp(ns_RESULT, 'ns_OK') ~= 1 )
    disp(['ERROR: ns_CloseFile() returned ' ns_RESULT]);
    return;
end

clear hFile
