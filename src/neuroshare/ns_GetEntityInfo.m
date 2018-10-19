function [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID)
% ns_GetEntityInfo - Retrieves general entity information and type
% 
% Usage:
% [ns_RESULT, nsEntityInfo] = ns_GetEntityInfo(hFile, EntityID)
% 
% Description:
% Retrieves general information about the entity, EntityID, from the file 
% referenced by the file handle hFile.  The information is passed in the 
% structure in nsEntityInfo. nsEntityInfo is a structure.
%
% Parameters:
% hFile                            A handle that contains information for 
%                                  one or more files. hFile is obtained by
%                                  a call to ns_OpenFile.  The result of 
%                                  ns_OpenFile, hFile, can be given as an 
%                                  argument in its entirety, or a subset of
%                                  the hFile array can be given.
%
% EntityID                         Identification number of the entity in 
%                                  the data file.
% 
% 
% Return Values:
% nsEntityInfo                     A cell array of ns_ENTITYINFO structures
%                                  to receive entity information for each 
%                                  of the files referenced in the hFile 
%                                  handle.
% 
%struct ns_ENTITYINFO
% char EntityLabel[32];          Specifies the label or name of the entity.
%
% double EntityType;             Flag specifying the type of entity data 
%                                recorded on this channel.
%                                It can be one of the following:
% 
%                                Unknown entity                   0
% 
%                                Event entity                     1
% 
%                                Analog entity                    2
% 
%                                Segment entity                   3
% 
%                                Neural event entity              4
% 
% double ItemCount;              Number of data items (samples) for the 
%                                specified entity in the file.
%
% ns_RESULT          This function returns one of the following status 
%                    codes:
%
%   ns_OK              The file was successfully opened. 
%   ns_BADFILE         Invalid file handle passed to function. 
%   ns_BADENTITY       Invalid or inappropriate entity identifier specified
%   ns_FILEERROR       File access or read error
%
% See also ns_OpenFile, ns_GetEntityInfo, ns_GetAnalogInfo,
% ns_GetAnalogData, ns_CloseFile

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
nsEntityInfo = [];

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
% list of entity types
EntityTypes = {'Event', 'Analog', 'Segment', 'Neural'};
ns_RESULT = 'ns_OK';

% start with empty label
nsEntityInfo.EntityLabel = [];

% if the Label field on the hFile Entity exists then copy it into the label on the nsEntityInfo 
if isfield(hFile.Entity(EntityID), 'Label')
    nsEntityInfo.EntityLabel = hFile.Entity(EntityID).Label;
    % IJM Added if logic below. 
    if size(nsEntityInfo.EntityLabel, 2) == 1
        nsEntityInfo.EntityLabel = [];
    end
end

% if the label on the nsEntityInfo is empty, auto-generate a label
if isempty(nsEntityInfo.EntityLabel)
    if strcmp(hFile.Entity(EntityID).EntityType, 'Analog')
        file_index = hFile.Entity(EntityID).FileType;
        sample_rate = 30/hFile.FileInfo(file_index).Period;
        nsEntityInfo.EntityLabel = sprintf('%d - %d kS/s',...
            hFile.Entity(EntityID).ElectrodeID, sample_rate);
    elseif strcmp(hFile.Entity(EntityID).EntityType, 'Event')
        % nsEntityInfo.EntityLabel = 'digin';
        nsEntityInfo.EntityLabel = hFile.Entity(EntityID).Reason;
    elseif ~isempty(hFile.Entity(EntityID).ElectrodeID)
        nsEntityInfo.EntityLabel =...
            ['elec' int2str(hFile.Entity(EntityID).ElectrodeID)];
    else
        nsEntityInfo.EntityLabel = hFile.Entity(EntityID).Reason;
    end
end

nsEntityInfo.EntityType = find(strcmp(EntityTypes,...
    hFile.Entity(EntityID).EntityType));
nsEntityInfo.ItemCount = hFile.Entity(EntityID).Count;
