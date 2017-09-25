function [ns_RESULT, entityIDs] = ...
    get_ripple_entity(hFile, port, headstage, channel)
% GET_RIPPLE_ENTITY - Retrieves analog data specific to analog entities 
% returns an array of structures the contain the Neuroshare index and 
% the entity type as a string for the wanted electrode.
% If no return value is specified writes wanted indices to stdout.
% Usage:
% [ns_RESULT, entityIDs] = get_ripple_entity(hFile, port, ...
%                                            headstage, channel)
% Parameters: 
%   hFile - Neuroshare file handle
%   port - letter of Grapevine port (a-d)
%   headstage - index of Grapevine headstage (1-4)
%   channel - index of front end pin for given headerstage (1-32)

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

% constants associated with the NIP
% Each lettered port chan hold a maximum of 128 channels
MAX_CHANNEL_PORT = 128; 
% each front-end can hold a max of 32 channels, up to 4 may be attached
% to a given lettered port
MAX_CHANNEL_HEADSTAGE = 32; 

entityIDs = [];

% ensure that we have a valid file handle
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return;
end

% allow for upper or lower case listed ports
port = lower(port);
% check that the port found is a single letter a-d
if strcmp(port, 'a')
    port_index = 0;
elseif strcmp(port, 'b')
    port_index = 1;
elseif strcmp(port, 'c')
    port_index = 2;
elseif strcmp(port, 'd')
    port_index = 3;
else
    fprintf(2, 'error: invalid NIP port %s\n', port);
    ns_RESULT = 'ns_BADENTITY';
    return;
end

% Map desired electrode to electrode_id.  This is valid only for 
% Ripple NIPs
elec_id = MAX_CHANNEL_PORT*port_index + ...
    (headstage-1)*MAX_CHANNEL_HEADSTAGE + channel;

% Find all the Neuroshare entities that have the desired electrode id
entities = find([hFile.Entity(:).ElectrodeID]==elec_id);
if isempty(entities)
    fprintf(2, 'error: could not find wanted channel\n');
    ns_RESULT = 'ns_BADENTITY';
    return
end

% Return an array of structures that contain the entity type,
% which is one of 'Neural', 'Analog', or 'Segment' and the index of
% the given entity in the Neuroshare API.
entityData = struct('type', '', 'index', 0);
% preallocate the returned array entityIDs
entityIDs = repmat(entityData, length(entities), 0);
% iterate through the found entities and save the data
for iEntity=1:length(entities)
    index = entities(iEntity);
    entityIDs(iEntity).type = hFile.Entity(index).EntityType;
    entityIDs(iEntity).index = index;
end
ns_RESULT = 'ns_OK';
% if no output is requested dump the entity ids to stdout
if nargout < 2
    fprintf('%7s %10s %14s\n', 'Index', 'Type', 'Label');
    for iEntity=1:length(entities)
        index = entities(iEntity);
        % file_index = hFile.Entity(index).FileType;
        type = hFile.Entity(index).EntityType;
        [ ns_RESULT, nsEntityInfo ] = ...
            ns_GetEntityInfo(hFile, index);
        label = nsEntityInfo.EntityLabel;
        fprintf(1, '%7d %10s %14s\n', index, type, label);
    end
end
