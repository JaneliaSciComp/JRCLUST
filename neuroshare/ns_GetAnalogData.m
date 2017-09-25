function [ns_RESULT, ContCount, Data] =...
    ns_GetAnalogData(hFile, EntityID, StartIndex, IndexCount, scaleFlag)
% ns_GetAnalogData - Retrieves analog data specific to analog entities 
% 
% Usage:
% [ns_RESULT, ContCount, Data] =...
%     ns_GetAnalogData(hFile, EntityID, StartIndex, IndexCount)
%
% Description:
% Returns the data values associated with the Analog Entity indexed
% EntityID in the file referenced by hFile. The index of the first data
% value is StartIndex and the requested number of data samples is given by
% IndexCount. The requested data values are returned in the variable Data.
% Although the samples in an analog entity are indexed, they are guaranteed
% to be continuous in time and may contain gaps between some of the indexes
% When the requested data is returned, ContCount contains the number of
% continuous data points present in the data (starting at StartIndex). if
% the index range specified by StartIndex and IndexCount contains invalid
% indexes, the function will return ns_BADINDEX. Data and ContCount are
% both a cell array of structures.
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
% StartIndex                       The index (sample number) to start 
%                                  analog data extraction. This field is 
%                                  ignored if the pointer is set to NULL.
%
% IndexCount                       Number of analog values to retrieve.
%
% Return Values:
% ContCount                        A cell array containing the number of 
%                                  continuous data values starting at 
%                                  StartIndex for each of the files 
%                                  referenced in the hFile handle.
%
% Data                             A cell array of double precision values 
%                                  with the analog data for each of the 
%                                  files referenced in the hFile handle.
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
ContCount = [];
Data = [];

% check input arguments
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

% check EntityID
if ~isnumeric(EntityID)||...
        (uint16(EntityID)~=EntityID)||...
        ~strcmp(hFile.Entity(EntityID).EntityType, 'Analog')
    ns_RESULT = 'ns_BADENTITY';
    return
end

% Check IndexCount
if ~isnumeric([StartIndex,IndexCount]) || IndexCount<1 || StartIndex<1
    ns_RESULT = 'ns_BADINDEX';
    return
end

ns_RESULT = 'ns_OK';

% create fileInfo structure
fileInfo = hFile.FileInfo(hFile.Entity(EntityID).FileType);

% calculate packet information
IndexTotal = min(StartIndex+IndexCount-1, hFile.Entity(EntityID).Count);
IndexCount = IndexTotal-StartIndex+1;
Data = zeros(IndexCount,1);
nPointAll = cumsum(fileInfo.TimeStamps(2,:));
StartPacket = find(StartIndex<=nPointAll, 1, 'first');
EndPacket = find(IndexTotal<=nPointAll, 1, 'first');
nPacket = length(StartPacket:EndPacket);
% create read size for each data Packet
if nPacket == 1
    PacketSize = IndexCount;
else
    PacketSize = fileInfo.TimeStamps(2,StartPacket:EndPacket);
    PacketSize(1) = PacketSize(1) - StartIndex + 1;
    PacketSize(end) = PacketSize(end) - (sum(PacketSize)-IndexCount);
end
% calculate offset
if strcmp(fileInfo.FileTypeID, 'NEUCDFLT')
    bytesPerPoint = 4;
else
    bytesPerPoint = 2;
end
bytesSkip = bytesPerPoint*length(fileInfo.ElectrodeList) - bytesPerPoint;
% Note: In the below, ismembc is a faster equivalent to find(list=val, 1).
% This is a bit obscure, but provides speed-up of 10x.
startDataLoc = StartPacket * 9 * (strcmp(fileInfo.FileTypeID, 'NEURALCD') || ...
    strcmp(fileInfo.FileTypeID, 'NEUCDFLT'));
%
% Mitch - ismembc2 does not work if the list is unsorted, which may happen
% if the analog input channels come first in the data file. find() has been
% exchanged for now both here and in ns_GetAnalogDataBlock for consistency
%
% offset = fileInfo.BytesHeaders +...
%     startDataLoc + ...
%     (bytesSkip+bytesPerPoint)*(StartIndex-1) +...
%     bytesPerPoint*ismembc2(hFile.Entity(EntityID).ElectrodeID,fileInfo.ElectrodeList) - bytesPerPoint;

offset = double(fileInfo.BytesHeaders) + ...
    double(startDataLoc) + ...
    double(bytesSkip + bytesPerPoint) * double(StartIndex-1) + ...
    double(bytesPerPoint) * find(fileInfo.ElectrodeList == hFile.Entity(EntityID).ElectrodeID,1,'first') - ...
    double(bytesPerPoint);

fseek(fileInfo.FileID, offset, -1);

% read data
CountList = [0, cumsum(PacketSize)];
for k = 1:nPacket
    if bytesPerPoint == 4
        Data(CountList(k)+1:CountList(k+1)) = ...
            fread(fileInfo.FileID, PacketSize(k), '*single', bytesSkip);
    else
        Data(CountList(k)+1:CountList(k+1)) = ...
            fread(fileInfo.FileID, PacketSize(k), '*int16', bytesSkip);
    end
    fseek(fileInfo.FileID, 9, 0);
end
% set size of returned data
ContCount = length(Data);
% if using float streams the data is already scaled
if bytesPerPoint == 4
    return;
end
% apply scale factor
% ContCount = PacketSize(1);

if exist('scaleFlag', 'var')
    if strcmpi(scaleFlag, 'unscale')
        return
    end
end
Data = Data*hFile.Entity(EntityID).Scale;
