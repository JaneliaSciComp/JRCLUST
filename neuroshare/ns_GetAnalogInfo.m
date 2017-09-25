function [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID)
% ns_GetAnalogInfo - Retrieves information specific to analog entities
% 
% Usage:
% [ns_RESULT, nsAnalogInfo] = ns_GetAnalogInfo(hFile, EntityID)
%
% Description:
% Returns information about the Analog Entity associated with EntityID and 
% the file hFile.  The information stored in nsAnalogInfo is a structure.
% Each index in the array refers to the entity
% information for each file in hFile.
%
% Parameters:
% hFile                       A handle that contains information for 
%                             one or more files. hFile is obtained by
%                             a call to ns_OpenFile.  The result of 
%                             ns_OpenFile, hFile, can be given as an 
%                             argument in its entirety, or a subset of
%                             the hFile array can be given.
% 
% EntityID                    Identification number of the entity in 
%                             the data file.
%
% Return Values:
% nsAnalogInfo               A cell array of ns_ANALOGINFO structures to 
%                            receive the Analog Entity information for each
%                            of the files referenced in the hFile handle.
% 
%struct ns_ANALOGINFO
% double SampleRate;          The sampling rate in Hz used to digitize the 
%                             analog values.
% 
% double MinVal;              Minimum possible value of the input signal.
% 
% double MaxVal;              Maximum possible value of the input signal.
% 
% char Units[16];             Specifies the recording units of measurement.
% 
% double Resolution;          Minimum input step size that can be resolved. 
% 
% double LocationX;           X coordinate of source in meters.
% 
% double LocationY;           Y coordinate of source in meters.
% 
% double LocationZ;           Z coordinate of source in meters.
% 
% double LocationUser;        Additional manufacturer-specific position 
%                             information.
% 
% double HighFreqCorner;      High frequency cutoff in Hz of the source 
%                             signal filtering.
% 
% double HighFreqOrder;       Order of the filter used for high frequency
%                             cutoff.
% 
% char HighFilterType[16];    Type of filter used for high frequency cutoff
%                             (text format).
% 
% double LowFreqCorner;       Low frequency cutoff in Hz of the source 
%                             signal filtering.
% 
% double LowFreqOrder;        Order of the filter used for low frequency 
%                             cutoff.
% 
% char LowFilterType[16];     Type of filter used for low frequency cutoff 
%                             (text format)..
% 
% char ProbeInfo[128];        Additional text information about the signal 
%                             source.
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
nsAnalogInfo = [];
% check hFile
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

%check Entity
if ~isnumeric(EntityID)||...
        (uint16(EntityID)~=EntityID)||...
        ~strcmp(hFile.Entity(EntityID).EntityType, 'Analog')
    ns_RESULT = 'ns_BADENTITY';
    return
end
ns_RESULT = 'ns_OK';
% create fileInfo structure for specified Entity
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);
% load default values of nsAnalogInfo structure
nsAnalogInfo.SampleRate = 30000/fileInfo.Period;
nsAnalogInfo.MinVal = [];
nsAnalogInfo.MaxVal = [];
nsAnalogInfo.Units = hFile.Entity(EntityID).Units;
nsAnalogInfo.Resolution = hFile.Entity(EntityID).Scale;
nsAnalogInfo.LocationX = [];
nsAnalogInfo.LocationY = [];
nsAnalogInfo.LocationZ = [];
nsAnalogInfo.LocationUser = [];
nsAnalogInfo.HighFreqCorner = [];
nsAnalogInfo.HighFreqOrder = [];
nsAnalogInfo.HighFilterType = [];
nsAnalogInfo.LowFreqCorner = [];
nsAnalogInfo.LowFreqOrder = [];
nsAnalogInfo.LowFilterType = [];
nsAnalogInfo.ProbeInfo = [];
% if nsx2.2 file load scale and filter information
if strcmp(fileInfo.FileTypeID, 'NEURALCD') || strcmp(fileInfo.FileTypeID, 'NEUCDFLT')
    FilterType = {'none', 'Butterworth', 'Chebyshev'};
    % find index of Entity Channel in Electrode List in order to obtain
    % specific extended header information
    % ELB_NOTE: I'm confused on this... Seeming to break with some files.
    % chanIDX = ismembc2(hFile.Entity(EntityID).ElectrodeID,...
    % fileInfo.ElectrodeList);
    chanIDX = find(fileInfo.ElectrodeList == hFile.Entity(EntityID).ElectrodeID);
    % calculate the start of the CC extended header
    headerPos = 314 + 66*(chanIDX - 1);
    % skip to ElectrodeLabel    
    fseek(fileInfo.FileID, headerPos + 4, -1);
    nsAnalogInfo.ProbeInfo = (fread(fileInfo.FileID, 16, '*char'))';
    % ELB, 2012/05/17
    % When the binary data is found in this label, The file pointer
    % reads one extra byte.  The rest of the read here is then garbage.
    % As a fix here, the pointer is advanced to where it's supposed to
    % from the absolute file position.
    fseek(fileInfo.FileID, headerPos + 20, -1);
    % skip: Physical Connector, Connector Pin, Min/Max Digital Value
    fseek(fileInfo.FileID, 6, 0);
    minMaxAnalog = fread(fileInfo.FileID, 2, 'int16');
    
    nsAnalogInfo.MinVal = minMaxAnalog(1);
    nsAnalogInfo.MaxVal = minMaxAnalog(2);
    % skip: Units
    fseek(fileInfo.FileID, 16, 0);
    highFilter = fread(fileInfo.FileID, 2, 'uint32');
    nsAnalogInfo.HighFreqCorner = highFilter(1);
    nsAnalogInfo.HighFreqOrder = highFilter(2);
    filterIDX = fread(fileInfo.FileID, 1, 'uint16');
    nsAnalogInfo.HighFilterType = char(FilterType(filterIDX+1));
    lowFilter = fread(fileInfo.FileID, 2, 'uint32');
    nsAnalogInfo.LowFreqCorner = lowFilter(1);
    nsAnalogInfo.LowFreqOrder = lowFilter(2);
    filterIDX = fread(fileInfo.FileID, 1, 'uint16');
    nsAnalogInfo.LowFilterType = char(FilterType(filterIDX+1));    
end
