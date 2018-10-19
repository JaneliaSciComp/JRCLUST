function [ns_RESULT, nsSegmentInfo] = ns_GetSegmentInfo(hFile, EntityID)
% ns_GetSegmentInfo - Retrieves information specific to segment entities
% Usage:
% 
% [ns_RESULT, nsSegmentInfo] = ns_GetSegmentInfo(hFile, EntityID)
% Description:
% 
% Retrieves information on the Segment Entity, EntityID, in the file 
% referenced by the handle hFile. The information is written to 
% nsSegmentInfo.
%
% Parameters:
% 
% hFile                      Handle/Identification number to an open file.
% 
% EntityID                   Identification number of the entity in the
%                            data file.
% Return Values:
% 
% nsSegmentInfo              ns_SEGMENTINFO structure that receives segment
%                            information for the requested Segment Entity.
% 
% struct ns_SEGMENTINFO
% 
% double SourceCount;        Number of sources contributing to the Segment 
%                            Entity data. For example, with tetrodes, this 
%                            number would be 4.
% 
% double MinSampleCount;     Minimum number of samples in each Segment data
%                            item.
% 
% double MaxSampleCount;     Maximum number of samples in each Segment data
%                            item.
% 
% double SampleRate;         The sampling rate in Hz used to digitize 
%                            source signals.
% 
% char Units[32];            Specifies the recording units of measurement.
% 
% ns_RESULT                  This function returns ns_OK if the file is 
%                            successfully opened. Otherwise one of the 
%                            following error codes is generated:
% 
% ns_BADFILE                 Invalid file handle passed to function
% 
% ns_BADENTITY               Invalid or inappropriate entity identifier 
%                            specified
%
% ns_FILEERROR               File access or read error

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
nsSegmentInfo = [];
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
ns_RESULT = 'ns_OK';
nsSegmentInfo.SourceCount = 1;
SampWaveforms =...
    (hFile.FileInfo(hFile.Entity(EntityID).FileType).BytesDataPacket-8)/2;
nsSegmentInfo.MinSampleCount = SampWaveforms;
nsSegmentInfo.MaxSampleCount = SampWaveforms;
nsSegmentInfo.SampleRate = 30000;
nsSegmentInfo.Units = hFile.Entity(EntityID).Units;
