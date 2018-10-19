% $Id$
% =========================================================================
%
%   PURPOSE: Offline data analysis tool set
% AUTHOR(S): Andrew Wilder
%   CONTACT: support@rppl.com
%
% (c) Copyright 2012 Ripple, LLC. All rights reserved.
%
% This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
% WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
%
% =========================================================================

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

%% Setup

close all;
clear all;

disp(' ');
disp(' ');
disp('Starting Trellis polling demo.');
disp(' ');

% misc variable initialization
sep = filesep;

% NEV extensions
nevExts = {'nev', 'ns1', 'ns2', 'ns3', 'ns4', 'ns5'};

% get the user's home directory
if ispc
    userHomeDir = getenv('USERPROFILE'); 
else
    userHomeDir= getenv('HOME');
end

% if this fails just start at the system root
if ~exist(userHomeDir, 'file')
    userHomeDir = '';
end

% select a directory to poll
pollDirPath = uigetdir(userHomeDir,'Select A Directory To Poll');

% see if the user has pressed "Cancel"
if pollDirPath == 0
    disp(' ');
    disp('Goodbye!');
    disp(' ');
    break;
end

% get the current contents of the poll dir
items = dir(pollDirPath);
processedDirExists = false;
for i=1:length(items)
    if (items(i).isdir) && (strcmp(items(i).name, 'processed'))
        processedDirExists = true;
        break;
    end
end

% if it does not exist create a directory for processed data files
if ~processedDirExists
    [success, message, messageID] = mkdir(pollDirPath, 'processed');
    if ~success
        disp(message);
        break;
    end
end
processedDirPath = [pollDirPath sep 'processed'];


%% Execute

disp('#');
disp('# Grapevine data processing script');
disp('#');
disp('#');

% data structure to keep trak of file size. need this to figure out whether 
% a NEV recording is complete
recordingSizes = containers.Map();

% poll the directory for new NEV recordings
while(1)
    
    % keep a list of recordings and their cumulative size (in bytes)
    latestSizes = containers.Map();
    
    % list the directory contents
    items = dir(pollDirPath);
    for i=1:length(items)
        
        % if the item is a file check if it is from a NEV recording
        if ~(items(i).isdir) 
            
            % parse the extension from the name
            [baseName, remain]  = strtok(items(i).name, '.');
            [ext,      remain2] = strtok(remain, '.');
            
            % check if this is a NEV recording file
            if sum(strcmp(ext, nevExts))
                
                % keep track of the file size
                recordingSize = items(i).bytes;

                % if we've seen other files from this recording compute
                % the sum of the sizes
                if( latestSizes.isKey(baseName) )
                    recordingSize = recordingSize + latestSizes(baseName);
                end
                    
                % keep track of the latest recording size
                latestSizes(baseName) = recordingSize;
            end
        end
    end
    
    % loop through the latest recording sizes and see if they match
    % previous recording sizes
    latestRecordings = latestSizes.keys;
    processedCount = 0;
    for i=1:latestSizes.Count
        
        % recording name
        recordingName = latestRecordings{i};
        
        % latest recording size
        latestSize = latestSizes(recordingName);
        
        % check if we've seen this recording before
        % and if the recording size is unchanged since the last check
        if     (recordingSizes.isKey(recordingName)) ...
           &&  (recordingSizes(recordingName) == latestSize) 
            
            % close all current plots
            close('all');
            
            % run the neuroshare demo
            disp(['Visualizing NEV recording: ' recordingName]);
            RippleNeuroshareDemo1(strcat(pollDirPath, sep, recordingName));
            
            processedCount = processedCount + 1;
                        
            % remove the recording from the list of recording sizes
            recordingSizes.remove(recordingName);
            
            % move all the files associated with this recording to 
            % the processed dir
            fromPath = strcat(pollDirPath, sep, recordingName, '*');
            toPath   = processedDirPath;
            [status, message, messageid] = movefile(fromPath, toPath);
            if ~status
                disp(['ERROR: could not move files from recording [' ...
                       recordingName ...
                      '] to processed dir.']);
            end
            
        % else keep track of the latest size of this recording
        else
            recordingSizes(recordingName) = latestSize;
        end
        
    end
    
    % generate status information
    dt       = fix(clock);
    dtStr    = [ num2str(dt(4), '%02d') ':' ...
                 num2str(dt(5), '%02d') ':' ...
                 num2str(dt(6), '%02d') ];
    countStr = num2str(processedCount);
    infoStr  = [ dtStr ' - No new NEV recordings'];
    
    % display status information
    fprintf(1,infoStr);
    
    % wait for a while before checking again for NEV recordings
    pause(1);
    
    % clear the status info
    for i=1:length(infoStr) fprintf(1,'\b'); end
end
