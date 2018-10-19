function plot_spikes(hFile, entityID, units)
%PLOT_SPIKES plot spike waveform to current pad.
%Usage: ns_RESULT = plot_spikes(hFile, entityID)
%   overlay all spike waveforms to current pad.
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       entityID - Wanted entity index, must be Segment entity
%       units - If specified, only plot spikes with unit ids in unit

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

% ensure that we have a valid NS file handle before we get started
if ~isstruct(hFile)
    fprintf(2, 'ns_BADFILE error\n');
    return
end
% Get ns_FileInfo and ns_EntityInfo structs.  If these functions
% fail return with error code
[ns_RESULT, fileInfo] = ns_GetFileInfo(hFile);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end
[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end
% Ensure that we have a Segment (spike) entity
if ~strcmp(hFile.Entity(entityID).EntityType,'Segment')
    fprintf(2, 'ns_BADENTITY error\n');
    return
end

% clear current subplot and prepare to overlay waveforms
reset(cla);
hold on;
% setup x-axis to have physical units.  Spike waveforms are
% always 52 bins.  x-axis is put in milliseconds

% Iterate through all items
for iItem=1:entityInfo.ItemCount
    [ns_RESULT, ts, data, sample_count, unit_id] = ...
        ns_GetSegmentData(hFile, entityID, iItem);
    % ts_res = fileinfo.TimeStampResolution;
    t = double(0:(sample_count-1))./30;
    
    % if there is something wierd about this entity, skip it
    if ~strcmp(ns_RESULT, 'ns_OK')
        fprintf(2, 'warning bad entity: %d\n', entityID);
        continue
    end
    % if wanted units are specified, only plot the units that
    % are contained in the array 'units'
    if nargin == 3
        if isscalar(units) && units ~= unit_id
            continue
        end
        if isempty(find(units==unit_id, 1))
            continue
        end
    end
    % color spike plots based on unit class.  These unit
    % classifications are compatible with Ripple data files
    color = 'black';
    if unit_id == 1
        color = 'red';
    elseif unit_id == 2
        color = 'blue';
    elseif unit_id == 3
        color = [0.23, 0.44, 0.34]; % dark green
    elseif unit_id == 4
        color = [0.68, 0.47, 0.0]; % deep yellow
    end
    % finally plot the waveforms
    plot(t, data, 'Color', color);
end
% setup axes and return to 'hold off'
title(entityInfo.EntityLabel);
xlabel('[ms]');
ylabel('[\muV]');
hold off;
