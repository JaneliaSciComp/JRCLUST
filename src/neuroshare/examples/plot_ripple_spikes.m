function plot_ripple_spikes(hFile, port, headstage, channel, varargin)
%PLOT_RIPPLE_SPIKE Plot spike data using Grapevine port ids
% This function is a wrapper for PLOT_SPIKES, but allows for parameters
% that translate easily to Ripple Grapevine NIP ports and headstages.
% 
% Usage: ns_RESULT = plot_ripple_spikes(hFile, port, headstage, channel)
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       port - single character that should be 'a' - 'd'
%       headstage - index of headstage used
%       channel - index of channel wanted relative to the headstage.
% See also PLOT_SPIKE

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

[ns_RESULT, entityIDs] = ...
    get_ripple_entity(hFile, port, headstage, channel);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end 
entity = entityIDs(strcmp({entityIDs(:).type}, 'Segment'));
index = entity.index;
if ~isempty(varargin)
    plot_spikes(hFile, index, varargin{1});
else
    plot_spikes(hFile, index);
end


