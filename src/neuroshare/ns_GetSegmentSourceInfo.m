function [ns_RESULT, nsSegmentSourceInfo] =...
    ns_GetSegmentSourceInfo(varargin)
% ns_GetSegmentSourceInfo - Retrieves information about the sources that 
%                           generated the segment data
% Usage:
% 
% [ns_RESULT, nsSegmentSourceInfo] =...
%             ns_GetSegmentSourceInfo(hFile, EntityID, SourceID)
% Description:
% 
% Retrieves information about the source entity, SourceID, for the Segment 
% Entity identified by EntityID, from the file referenced by the handle 
% hFile. The information is written to the nsSegmentSourceInfo.
%
% Parameters:
% 
% hFile                     Handle/Identification number to an open file.
% 
% EntityID                  Identification number of the Segment Entity.
% 
% SourceID                  Identification number of the Segment Entity 
%                           source (integer from 0 to SourceCount – 1 which
%                           can be found in ns_SEGMENTINFO structure).
%
% Return Values:
% 
% nsSegmentSourceInfo       ns_SEGSOURCEINFO structure that receives 
%                           information about the source.
% 
% struct ns_SEGSOURCEINFO
% 
% double MinVal;            Minimum possible value of the input signal.
% 
% double MaxVal;            Maximum possible value of the input signal.
% 
% double Resolution;        Minimum input step size that can be resolved. 
% 
% double SubSampleShift;    Time difference (in sec) between the nominal
%                           timestamp and the actual sampling time of the
%                           source probe.
% 
% double LocationX;         X coordinate of source in meters.
% 
% double LocationY;         Y coordinate of source in meters.
% 
% double LocationZ;         Z coordinate of source in meters.
% 
% double LocationUser;      Additional manufacturer-specific position 
%                           information
% 
% double HighFreqCorner;    High frequency cutoff in Hz of the source 
%                           signal filtering.
% 
% double HighFreqOrder;     Order of the filter used for high frequency 
%                           cutoff.
% 
% char HighFilterType[16];  Type of filter used for high frequency cutoff 
%                           (text format).
% 
% double LowFreqCorner;     Low frequency cutoff in Hz of the source signal
%                           filtering.
% 
% double LowFreqOrder;      Order of the filter used for low frequency 
%                           cutoff.
% 
% char LowFilterType[16];   Type of filter used for low frequency cutoff 
%                           (text format).
% 
% char ProbeInfo[128];      Additional text information about the signal 
%                           source.
% 
% ns_RESULT                 This function returns ns_OK if the file is 
%                           successfully opened. Otherwise one of the 
%                           following error codes is generated:
% 
% ns_BADFILE                Invalid file handle passed to function
% 
% ns_BADENTITY              Invalid or inappropriate entity identifier 
%                           specified
% 
% ns_BADSOURCE              Invalid source identifier specified
% 
% ns_FILEERROR              File access or read error

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

narginchk(2, 3);
hFile = varargin{1};
EntityID = varargin{2};

nsSegmentSourceInfo = [];
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

if nargin == 2
    SourceID = 1;
else
    SourceID = varargin{3};
end

%check SourceID, must be 1 for nev/nsx files
if SourceID ~=1
    ns_RESULT = 'ns_BADSOURCE';
    return
end

ns_RESULT = 'ns_OK';

nsSegmentSourceInfo.MinVal = [];
nsSegmentSourceInfo.MaxVal = [];
nsSegmentSourceInfo.Resolution = hFile.Entity(EntityID).Scale;
nsSegmentSourceInfo.SubSampleShift = [];
nsSegmentSourceInfo.LocationX = [];
nsSegmentSourceInfo.LocationY = [];
nsSegmentSourceInfo.LocationZ = [];
nsSegmentSourceInfo.LocationUser = [];
nsSegmentSourceInfo.HighFreqCorner = [];
nsSegmentSourceInfo.HighFreqOrder = [];
nsSegmentSourceInfo.HighFilterType = [];
nsSegmentSourceInfo.LowFreqCorner = [];
nsSegmentSourceInfo.LowFreqOrder = [];
nsSegmentSourceInfo.LowFilterType = [];
nsSegmentSourceInfo.ProbeInfo = ...
    int2str(hFile.Entity(EntityID).ElectrodeID);

wantedElectrode = hFile.Entity(EntityID).ElectrodeID;
% get file handle for this entities nev file
fid = hFile.FileInfo(hFile.Entity(EntityID).FileType).FileID;
% Find the NEUEVFLT header that has the info for this electrode
% NOTE: This method was  more or less copied this method from ns_OpenFile. 
% Should be improved if it proves to be slow.
fseek(fid, 332, -1);
nExtendedHeaders = fread(fid, 1, '*uint32');

% Get list of extended header PacketIDs: (there are 9 different nev
% extended headers. Each have a 8 char array PacketID and each
% extended header has 24 bytes of information with a total size of
% 8+24=32 bytes. 
PacketIDs = cellstr(fread(fid, [8, nExtendedHeaders],...
              '8*char=>char', 24)');
          
% set filter types
FilterType = {'none', 'Butterworth', 'Chebyshev'};
          
% Get Index of NEUEVWAV extended headers.
idxEVFLT = find(strcmp(PacketIDs, 'NEUEVFLT'));

for j = 1:length(idxEVFLT)
    % seek to the next NEUEVWAV extended header from the begining
    % of the file. 
    fseek(fid, 344+(idxEVFLT(j)-1)*32, -1);
    elecID = fread(fid, 1, '*uint16');
    %fprintf('%d\n', elecID);
    if elecID == wantedElectrode
        nsSegmentSourceInfo.HighFreqCorner = fread(fid, 1, 'uint32');
        nsSegmentSourceInfo.HighFreqOrder = fread(fid, 1, '*uint32');
        highType = fread(fid, 1, '*uint16');
        nsSegmentSourceInfo.HighFilterType = char(FilterType(highType+1));
        nsSegmentSourceInfo.LowFreqCorner = fread(fid, 1, 'uint32');
        nsSegmentSourceInfo.LowFreqOrder = fread(fid, 1, '*uint32');        
        lowType = fread(fid, 1, '*uint16');
        nsSegmentSourceInfo.LowFilterType = char(FilterType(lowType+1));
        ns_RESULT = 'ns_OK';
        return
    end
end

ns_RESULT = 'ns_BADENTITY';
return
