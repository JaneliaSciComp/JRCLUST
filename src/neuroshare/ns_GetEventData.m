function [ns_RESULT, TimeStamp, Data, DataSize] = ...
    ns_GetEventData(hFile, EntityID, Index)
% ns_GetEventData - Retrieves event data by index
% Usage:
% [ns_RESULT, TimeStamp, Data, DataSize] =...
%                        ns_GetEventData(hFile, EntityID, Index)
% Description:
% 
% Returns the data values from the file referenced by hFile and the Event 
% Entity EntityID.  The Event data entry specified by Index is written to 
% Data and the timestamp of the entry is returned to TimeStamp.  
% Upon return of the function, the value at DataSize contains the number of
% bytes actually written to Data.
%
% Parameters:
%
% hFile               Handle/Identification number to an open file.
% EntityID            Identification number of the entity in the data file.
% Index               The index number of the requested Event data item.
%
% Return Values:
% 
% TimeStamp           Variable that receives the timestamp of the Event 
%                     data item.
%
% Data                Variable that receives the data for the Event entry.
%                     The format of data is specified by the member 
%                     EventType in ns_EVENTINFO.
%
% DataSize            Variable that receives the actual number of bytes of 
%                     data retrieved in the data buffer.
%
% ns_RESULT           This function returns ns_OK if the file is 
%                     successfully opened. Otherwise one of the following 
%                     error codes is generated:
%
% ns_BADFILE          Invalid file handle passed to function
% ns_BADENTITY        Invalid or inappropriate entity identifier specified
% ns_BADINDEX         Invalid entity index specified
% ns_FILEERROR        File access or read error

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

TimeStamp = [];
Data = [];
DataSize = [];
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

%check Entity
if ~(isnumeric(EntityID))||...
    (uint16(EntityID)~=EntityID)||...
    (EntityID > length(hFile.Entity))||...
    (EntityID < 1)||...
    ~strcmp(hFile.Entity(EntityID).EntityType, 'Event')
    ns_RESULT = 'ns_BADENTITY';
    return
end

% check Index
if Index<1||Index>hFile.Entity(EntityID).Count
    ns_RESULT = 'ns_BADINDEX';
    return
end

ns_RESULT = 'ns_OK';
DataSize = 2;
% construct the fileInfo structure for the specifies entityID
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);
% Use this cell array to match the reason for this digital event
% to be stored.  This string is held internally
% packetReason = {'Digital Input',...
%                         'Input Ch 1',...
%                         'Input Ch 2',...
%                         'Input Ch 3',...
%                         'Input Ch 4',...
%                         'Input Ch 5'};
% IJM - Modified to support new labeling scheme. 
% Should I rewrite the code to support the legacy labels? 
packetReason = {'Parallel Input', 'SMA 1', 'SMA 2', 'SMA 3', ...
                    'SMA 4', 'Output Echo', 'Digital Input', ...
                    'Input Ch 1', 'Input Ch 2', 'Input Ch 3', ...
                    'Input Ch 4', 'Input Ch 5'};
% Using the cell array find the get the bit at determines the
% channel crossing for this entity to be stored
idx = find(strcmp(packetReason,  hFile.Entity(EntityID).Reason));
% IJM - Adjust index for the old packet reasons.
if idx > 6
    idx = idx - 6;
end
% Create a list of all the digital events with the classification for
% this digital entity.  This means, packetid==0 and  the stored Class
% has an up bit in position idx.
posEvent = ...
    find(bitget(fileInfo.MemoryMap.Data.Class, idx) == 1 & ...
        fileInfo.MemoryMap.Data.PacketID==0, Index);
% The because we ask only for Index values, the end of the above array
% is our desired entity

% Get the timestamp for this event.  Use a constant 30kHz sampling rate 
% which should be consitent for all nev files
TimeStamp = double(fileInfo.MemoryMap.Data.TimeStamp(posEvent(end)))/30000;

% calculate the offset in bytes to the data that we are interested in 
% from the start of the nev file
offset = double(fileInfo.BytesHeaders) +...
    double(fileInfo.BytesDataPacket) *...
    double(posEvent(end)-1) +...
    8 + (idx-1)*2;
% skip to the calculated offset
fseek(fileInfo.FileID, offset, -1);
% if the wanted data is the digital parrelle port only 
% unsigned ints are supported, otherwise expect signed ints
if idx == 1
    Data = fread(fileInfo.FileID, 1, '*uint16');
else
    Data = fread(fileInfo.FileID, 1, '*int16');
end
