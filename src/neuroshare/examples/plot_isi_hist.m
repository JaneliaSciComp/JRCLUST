function ns_RESULT = plot_isi_hist(varargin)
%PLOT_ISI_HIST Create ISI Histogram for this entity.
%Usage: ns_RESULT = plot_spikes(hFile, entityID, [bins|binSize])
%   If no binning is provided then the bin size used is 1 ms with
%   histogram boundries of 0 to 100 ms.
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       entityID - Wanted entity index, must be Segment entity
%   Optional:
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
% Require two parameters to specify NS file handle and wanted entity.
% narginchk is only valid in matlab versions later than 2012
% narginchk(2, 4);
if nargin < 2 || nargin > 4
    fprintf(2, 'error: invalid number of parameters.\n');
    return
end
hFile = varargin{1};
entityID = varargin{2};

% set up bins array
if nargin > 2
    if length(varargin{3}) > 1
        % if the third parameter is specified as an array,
        % use this as the bin centers.
        bins = varargin{3};
    else
        % if the third parameter is a scalar, use this 
        % as the bin centers.
        bins = 0:varargin{3}:100;
    end
else
    % if no bin size is specified use 1 ms as bin size and
    % make the histogram range 0 - 100ms
    bins = 0:1:100;
end

% Get ns_EntityInfo struct.  If this fails, return exit code
[ns_RESULT, entityInfo] = ns_GetEntityInfo(hFile, entityID);
if ~strcmp(ns_RESULT, 'ns_OK')
    return
end

% return a failure if we don't have a segment entity
if ~strcmp(hFile.Entity(entityID),'Segment') == 0
    ns_RESULT = 'ns_BADENTITY';
    return
end

% clear current subplot and prepare to draw hist
reset(cla);

% set up array to hold the inter-site interval values
isi = zeros(entityInfo.ItemCount-1);

% save the previous value from the loop
tsLast = -1;
% Iterate through all items and fill wanted timestamps
for iItem=1:entityInfo.ItemCount
    [ns_RESULT, ts, data, sample_count, unit_id] = ...
        ns_GetSegmentData(hFile, entityID, iItem);
    % if this fails, something wierd is happening, lets stop.
    if ~strcmp(ns_RESULT, 'ns_OK')
        return
    end
    % for the first value, store timestamp and continue.
    if tsLast == -1
        tsLast = ts;       
        continue
    end
    % set difference between this iteration and the last in msec
    isi(iItem-1) = (ts-tsLast)*1000.0;
    tsLast = ts;
end
isi = isi(isi~=0.0);
hist(isi, bins);
    title(entityInfo.EntityLabel);
xlim([min(bins) max(bins)]);
xlabel('[ms]');

ns_RESULT = 'ns_OK';
