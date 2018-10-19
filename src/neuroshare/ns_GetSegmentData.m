function [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ...
    ns_GetSegmentData(hFile, EntityID, Index)
% ns_GetSegmentData - Retrieves segment data by index
% Usage:
% 
% [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ...
%             ns_GetSegmentData(hFile, EntityID, Index)
% Description:
% 
% Returns the Segment data values in entry Index of the entity EntityID 
% from the file referenced by hFile. The data values are returned in Data. 
% The timestamp of the entry id returned in TimeStamp.  The number of 
% samples written to Data is returned in SampleCount.  The data variable 
% should be accessed as a 2-dimensional array for samples and sources.
% 
% Parameters:
% 
% hFile               Handle/Identification number to an open file.
% 
% EntityID            Identification number of the entity in the data file.
% 
% Index               Index number of the requested Segment data item.
%       
% Remarks:
% 
% A zero unit ID is unclassified, then follow unit 1, 2, 3, etc.
%
% Return Values:
% 
% TimeStamp           Time stamp of the requested Segment data item.
% 
% Data                Variable to receive the requested data.
% 
% SampleCount         Number of samples returned in the data variable.
% 
% UnitID              Unit classification code for the Segment Entity.
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

TimeStamp    = [];
Data         = [];
SampleCount  = [];
UnitID       = [];

%check hFile
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

%check Entity
if ~(isnumeric(EntityID))||...
    (uint16(EntityID)~=EntityID)||...
    (EntityID > length(hFile.Entity))||...
    (EntityID < 1)||...
    ~strcmp(hFile.Entity(EntityID).EntityType, 'Segment')
    ns_RESULT = 'ns_BADENTITY';
    return
end

%check Index
if Index<1||Index>hFile.Entity(EntityID).Count
    ns_RESULT = 'ns_BADINDEX';
    return
end

ns_RESULT = 'ns_OK';

%create fileInfo Structure
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);
SampleCount = (fileInfo.BytesDataPacket - 8)/2;
PacketIndex = find(fileInfo.MemoryMap.Data.PacketID == ...
    hFile.Entity(EntityID).ElectrodeID, Index);
TimeStamp = double(fileInfo.MemoryMap.Data.TimeStamp(PacketIndex(end)))/30000;
UnitID = fileInfo.MemoryMap.Data.Class(PacketIndex(end));
offset = double(fileInfo.BytesHeaders) + 8 + ...
    double(fileInfo.BytesDataPacket) * double(PacketIndex(end)-1);
            
% skip to event wave data
fseek(fileInfo.FileID, offset, -1);
Data = fread(fileInfo.FileID, SampleCount, 'int16=>double');

% scale the data
Data = Data*hFile.Entity(EntityID).Scale;
