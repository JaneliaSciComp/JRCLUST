function [ns_RESULT, nsEventInfo] = ns_GetEventInfo(hFile, EntityID)
% ns_GetEventInfo - Retrieves information specific to event entities
% Usage:
% 
% [ns_RESULT, nsEventInfo] = ns_GetEventInfo(hFile, EntityID)
% Description:
% 
% Retrieves information from the file referenced by hFile about the Event 
% Entity, EntityID, in the structure nsEventInfo.
%
% Parameters:
% struct nsEventInfo
%
% hFile                       Handle/Identification number to an open file.
% EntityID                    Identification number of the entity in the 
%                             data file.
%
% Return Values:
% 
% nsEventInfo                 ns_EVENTINFO structure to receive the Event 
%                             Entity information.
%
% double EventType;           A type code describing the type of event data 
%                             associated with each indexed entry.  
%                             The following information types are allowed:
% 
%                             Text string                   0
%                             Comma separated values        1
%                             8-bit binary values           2
%                             16-bit binary values          3
%                             32-bit binary values          4
% 
% double MinDataLength;       Minimum number of bytes that can be returned 
%                             for an Event.
% 
% double MaxDataLength;       Maximum number of bytes that can be returned 
%                             for an Event.
% 
% ns_RESULT                   This function returns ns_OK if the file is 
%                             successfully opened. Otherwise one of the 
%                             following error codes is generated:
% 
% ns_BADFILE                  Invalid file handle passed to function
% 
% ns_BADENTITY                Invalid or inappropriate entity identifier 
%                             specified
% 
% ns_FILEERROR                File access or read error

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
nsEventInfo = [];
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
    ~strcmp(hFile.Entity(EntityID).EntityType, 'Event')
    ns_RESULT = 'ns_BADENTITY';
    return
end
ns_RESULT = 'ns_OK';
nsEventInfo.EventType = 3;
nsEventInfo.MinDataLength = 2;
nsEventInfo.MaxDataLength = 2;
nsEventInfo.CSVDesc = [];