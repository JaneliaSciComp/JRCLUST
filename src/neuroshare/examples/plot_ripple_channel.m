function plot_ripple_channel(hFile, port, headstage, channel)
% PLOT_RIPPLE_CHANNEL Produce three panel plot for Ripple channel
% This function provides an easy interface to plot spike waveforms,
% spike times, and analog wave forms for the given Ripple channel
%
% Usage: ns_RESULT = plot_ripple_channel(hFile, port, headstage, channel)
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       port - single character that should be 'a' - 'd'
%       headstage - index of headstage used
%       channel - index of channel wanted relative to the headstage.
% See also PLOT_RIPPLE_SPIKE, PLOT_RIPPLE_HIST, PLOT_RIPPLE_ANALOG

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

% first just check that we have a valid input arguments by calling
% get_ripple_entity
[ns_RESULT, entityIDs] = get_ripple_entity(hFile, port,...
                                           headstage, channel);
if ~strcmp(ns_RESULT, 'ns_OK')
    % just exit, get_ripple_entity should already have outputted an error.
    return
end

% Get screen dimensions so that the sizes of the figure windows can
% be controlled for any screen resolution
SCREEN_SIZE = get(0, 'ScreenSize');
% modifiying header and footer parameters may beautify window sizes
SCREEN_FOOTER = 55;
SCREEN_HEADER = 150;
SCREEN_L_MARGIN = 5;

% Pack the wanted positions into an array.  This figure will take up
% the left 1/3 of the width of the screen
position = [ SCREEN_L_MARGIN SCREEN_FOOTER ...
    SCREEN_SIZE(3)/3 0.75*SCREEN_SIZE(4) ];

% if this function has been run before and the figure exists
% close the figure and start fresh
if exist('rippleHandle', 'var')
    close(rippleHandle)
end
rippleHandle = figure('Position', position, 'NumberTitle', 'off', ...
    'Name', 'Ripple Plot');

% plot spikes
subplot(3, 1, 1);
plot_ripple_spikes(hFile, port, headstage, channel);
% plot timestamp hist
subplot(3, 1, 2);
plot_ripple_hist(hFile, port, headstage, channel);
% plot analog data
subplot(3, 1, 3);
plot_ripple_analog(hFile, port, headstage, channel);
