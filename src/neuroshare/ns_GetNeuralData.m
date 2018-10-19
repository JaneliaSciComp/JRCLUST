function [ns_RESULT, Data] = ...
    ns_GetNeuralData(varargin)
%NS_GETNEURALDATA - Retrieves neural event timestamps in seconds
% Usage:
% 
% [ns_RESULT, TimeStamp, Data, SampleCount, UnitID] = ...
%             ns_GetSegmentData(hFile, EntityID, StartIndex, IndexCount)
% Description:
% 
% Returns the timestamps in seconds for Neural entities corresponding
% to the entity index EntityID.  Neural entity will be created for each
% trigger spike and each sorted event.  This functions allows for the 
% ability to get trigger times based on spike classifications.
% 
% Parameters:
% 
% hFile               Handle/Identification number to an open file.
% 
% EntityID            Identification number of the entity in the data file.
% 
% StartIndex          (Optional) Return timestamps after StartIndex 
%                     occurences.
%
% IndexCount          (Optional) Return only IndexCount timestamps
%                     


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

Data = [];

narginchk(2, 4);
%nargchk(2, 4);
hFile = varargin{1};
EntityID = varargin{2};

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
    ~strcmp(hFile.Entity(EntityID).EntityType, 'Neural')
    ns_RESULT = 'ns_BADENTITY';
    return
end
% The neural item count will be useful for checking the validity
% of StartIndex and IndexCount
itemCount = hFile.Entity(EntityID).Count;
if nargin >= 3
    StartIndex = varargin{3};
    if StartIndex < 1 || StartIndex > hFile.Entity(EntityID).Count
        ns_RESULT = 'ns_BADINDEX';
        return
    end
else
    % if StartIndex is not specified start at the beginning of the data
    StartIndex = 1;
end
if nargin == 4
   IndexCount = varargin{4}; 
   if IndexCount > itemCount
       ns_RESULT = 'ns_BADINDEX';
       return
   end
   endBin = StartIndex + IndexCount - 1;
else
   % if IndexCount is not specified go the end of the neural entities
   endBin = itemCount;
end
ns_RESULT = 'ns_OK';

% get needed variables out of the hFile struct
elecID = hFile.Entity(EntityID).ElectrodeID;
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);
class = hFile.Entity(EntityID).Reason;
% find instances where we have the wanted electrode id and sorted
% neural entity
neuralIndices = logical(fileInfo.MemoryMap.Data.PacketID == elecID &...
         fileInfo.MemoryMap.Data.Class == class);
% Get all the timestamps that match our criteria and put them in seconds
timestamps = ...
    double(fileInfo.MemoryMap.Data.TimeStamp(neuralIndices))/30000;
% slice the wanted part of the timestamps
Data = timestamps(StartIndex:endBin);
end
