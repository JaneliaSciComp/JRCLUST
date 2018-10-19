function [ns_RESULT, Time] = ns_GetTimeByIndex(hFile, EntityID, Index)
% ns_GetTimeByIndex - Retrieves time range from entity indexes
% Usage:
% 
% [ns_RESULT, Time] = ns_GetTimeByIndex(hFile, EntityID, Index)
% Description:
% 
% Retrieves the timestamp for the entity identified by EntityID 
% and numbered Index, from the data file referenced by hFile. The timestamp
% is returned in Time.
% 
% Parameters:
% 
% hFile               Handle/Identification number to an open file.
% 
% EntityID            Identification number of the entity in the data file.
% 
% Index               Index of the requested data.
% 
% Return Values:
% 
% Time                Variable to receive the timestamp.
% 
% ns_RESULT           This function returns ns_OK if the file is 
%                     successfully opened. Otherwise one of the following 
%                     error codes is generated:
% 
% ns_BADFILE          Invalid file handle passed to function
% 
% ns_BADENTITY        Invalid or inappropriate entity identifier specified
% 
% ns_BADINDEX         Invalid entity index specified
% 
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
Time = [];

%check hFile
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

%check Entity
if ~(isnumeric(EntityID))||...
    (EntityID > length(hFile.Entity))||...
    (EntityID < 1)||...
    (uint16(EntityID)~=EntityID)
    ns_RESULT = 'ns_BADENTITY';
    return
end

if ~(isnumeric(Index))||...
        Index<1||...
        Index > hFile.Entity(EntityID).Count
    ns_RESULT = 'ns_BADINDEX';
    return
end
ns_RESULT = 'ns_OK';
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);

if strcmp(hFile.Entity(EntityID).EntityType, 'Analog')
    nPointsAll = cumsum(fileInfo.TimeStamps(2,:));
    idx = find(Index<=nPointsAll,1);
    IndexList = [0, nPointsAll(1:end-1)]+1;
    Time = (fileInfo.TimeStamps(1,idx)+Index-IndexList(idx))...
        *fileInfo.Period/30000;
    return
elseif strcmp(hFile.Entity(EntityID).EntityType, 'Segment')
    TimeStamps = fileInfo.MemoryMap.Data.TimeStamp(...
        fileInfo.MemoryMap.Data.PacketID==...
        hFile.Entity(EntityID).ElectrodeID);
    Time = double(TimeStamps(Index))/30000;
    return
elseif strcmp(hFile.Entity(EntityID).EntityType, 'Event')
    packetReason = {'Digital Input',...
                            'Input Ch 1',...
                            'Input Ch 2',...
                            'Input Ch 3',...
                            'Input Ch 4',...
                            'Input Ch 5'};
    idx = find(strcmp(packetReason,  hFile.Entity(EntityID).Reason));
    %eventClass = fileInfo.MemoryMap.Data.Class(...
    %    fileInfo.MemoryMap.Data.PacketID == 0);
    %TimeStamps = fileInfo.MemoryMap.Data.TimeStamp(bitget(eventClass,idx));
    %Time = double(TimeStamps(Index))/30000;
    % IJM: Removed the above which gave a fault result. Took the snippet
    % below from ns_GetEventData. 
    posEvent = ...
    find(bitget(fileInfo.MemoryMap.Data.Class, idx) == 1 & ...
        fileInfo.MemoryMap.Data.PacketID==0, Index);
    Time = double(fileInfo.MemoryMap.Data.TimeStamp(posEvent(end)))/30000;
end
    
    