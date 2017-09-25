function varargin = plot_analog(varargin)
%PLOT_ANALOG plot Neuroshare analog data.
% Usage:
% ns_RESULT = plot_analog(hFile, entityID, [startIndex, indexCount])
%   plots the analog data for Neuroshare file handle hFile and entity
%   index entityID to the current global pad.
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       entityID - Wanted entity index, must be Analog entity
%   Optional: 
%       startIndex: first wanted bins (default 0)
%       indexCount: number of wanted bins to plot (default all)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The Wisteria Neuroshare Importer is free software: you can 
%     redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     The Wisteria Neuroshare Importer is distributed in the hope that it 
%     will be useful, but WITHOUT ANY WARRANTY; without even the implied 
%     warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
%     See the GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with the Wisteria Neuroshare Importer.  If not, see 
%     <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% handle a variable number of inputs
if nargin < 2 || nargin > 4
    fprintf(2, 'error: invalid number of parameters.\n');
    return
end

% hFile and entityID are required for this script to run
hFile = varargin{1};
entityID = varargin{2};

% startIndex and indexCount are optional.  These are the defaults
startIndex = 1;
indexCount = -1;
% if specified set startIndex and indexCount
if nargin > 2
    startIndex = varargin{3};
end
if nargin == 4
    indexCount = varargin{4};
end  

% ensure that we have a valid NS file handle before we get started
if ~isstruct(hFile)
    fprintf(2, 'ns_BADFILE error\n');
    return
end
% get ns_AnalogInfo and ns_EntityInfo structs.  If these functions
% return an error, exit this function with error code.
[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end
[ns_RESULT, analogInfo] = ns_GetAnalogInfo(hFile, entityID);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end
% ensure that we are working with an analog entity
if entityInfo.EntityType ~= 2
    fprintf(2, 'ns_BADENTITY error\n');
    return
end
% if no index count is specified use the total number of points 
if indexCount < 0
    indexCount = entityInfo.ItemCount;
end

% clear current subplot and prepare to draw hist
reset(cla);
% get analog data
[ns_RESULT, countList, data] = ...
    ns_GetAnalogData(hFile, entityID, startIndex, indexCount);

% count list contains and array of counts for each recording section
% if there are pauses.  Here we just add them up and put them on one plot.
count = 0;
for i=1:length(countList)
    count = count + countList(1);
end

% create physical time in seconds for plotting
t = (startIndex-1:startIndex+count-2)*1.0/analogInfo.SampleRate;

% exit with error code if we failed to find data
if ~strcmp(ns_RESULT, 'ns_OK')
    return
end
% finally plot the data
plot(t, data);

% setup titles and axis
title(entityInfo.EntityLabel);
xlabel('[s]');
ylabel('[\muV]');
xlim([0, max(t)]);
