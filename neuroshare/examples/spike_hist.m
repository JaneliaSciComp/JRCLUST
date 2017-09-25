function [ vargargout ] = spike_hist(varargin)
%SPIKE_HIST histogram spike times.
% Usage:
% [h, bins] = plot_analog(hFile, entityID, [units], [bins|binSize])
%   Plot or return histogram data for spike data.  Behaviour is identical
%   to hist.  
%   Returns:
%       nargout == 0 - return nothing, plot histogram
%       nargout == 1 - return bin contents
%       nargout == 2 - return bin contents and bin centers 
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       entityID - Wanted entity index, must be Segment entity
%   Optional:
%       units - If specified, only return these unit ids
%       bins or binSize - If this paramter is a scalar then it will be
%           taken as the bin size and the default range of 0-100ms is used.
%           If the the third parameter is an array, then it is used to
%           provide the bin centers for the histogram.

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

% handle a variable number of inputs
% narginchk is only valid in matlab versions later than 2012
% narginchk(2, 4);
if nargin < 2 || nargin > 4
    fprintf(2, 'error: invalid number of parameters.\n');
    return
end
% hFile and entityID are required for this script to run
hFile = varargin{1};
entityID = varargin{2};

if ~isstruct(hFile)
    fprintf(2, 'ns_BADFILE error\n');
    return
end
if ~strcmp(hFile.Entity(entityID).EntityType,'Segment')
    fprintf(2, 'ns_BADENTITY error\n');
    return
end
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

if nargin > 2
    if length(varargin{3}) > 1
        bins = varargin{3};
        binSize = 0;
    else
        binSize = varargin{3};
    end
else
    binSize = 0.2; % 200 millisec
end

% clear current subplot and prepare to draw hist
reset(cla);

tsData = struct('ts', 0.0, 'unit', -1);
timestamps = repmat(tsData, entityInfo.ItemCount, 1);
count = 0;
% Iterate through all items and fill wanted timestamps

for iItem=1:entityInfo.ItemCount
    [ns_RESULT, ts, data, sample_count, unit_id] = ...
        ns_GetSegmentData(hFile, entityID, iItem);
    % if wanted units are specified, only plot the units that
    % are contained in the array 'units'
    if nargin == 4
        units = varargin{4};
        if isscalar(units) && units~=unit_id
            continue
        end
        if isempty(find(units==unit_id, 1))
            continue
        end
    end
    count = count+1;
    timestamps(count).ts = ts;
    timestamps(count).unit = unit_id;
end
max_time = max([timestamps(:).ts]);
if binSize > 0
    bins = 0:binSize:max_time;
end
timestamps = timestamps([timestamps(:).unit]~=-1);
timestamps = timestamps([timestamps(:).ts]~=0);

% if output is specified return fill the histogram 
% and return all the wanted 
if nargout > 0
    vargargout = hist([timestamps(:).ts], bins);
    return
else
    hold on;
    % is_plot decides whether we are are going to plot or not
    is_plot = 0;
    for iSort=0:4
        good_times = timestamps([timestamps(:).unit] == iSort);
        ts = [good_times(:).ts];
        if isempty(ts)
            continue
        end
        [y, x] = hist(ts, bins);
        % color spike plots based on unit class.  These unit
        % classifications are compatible with Ripple data files
        color = 'black';
        if iSort == 1
            color = 'red';
        elseif iSort == 2
            color = 'blue';
        elseif iSort == 3
            color = [0.23, 0.44, 0.34]; % dark green
        elseif iSort == 4
            color = [0.68, 0.47, 0.0]; % deep yellow
        end
        bar(x, y, 'EdgeColor', color, 'FaceColor', color);
        is_plot = 1;
    end
end
hold off;
if is_plot
    % set(gca, 'ytick', []);
    % ensure that the x-axis times match the 
    % timestamps that were found in the entities.
    xlim([0, max_time]);
    title(entityInfo.EntityLabel);
    xlabel('[s]');
else
    fprintf('warning: did not find any valid spike times\n')
end
