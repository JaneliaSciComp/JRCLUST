function [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile)
% ns_GetFileInfo - Retrieves file information and entity counts
% 
% Usage:
% [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile)
% Returns a cell array containing data structures with file information for
% files referenced in the hFile argument.
%
% Description: 
% Provides general information about the data file(s) referenced by hFile. 
% ns_GetFileInfo returns nsFileInfo, a structure associated with the
% open files.
%
% Parameters:
% hFile                            A handle that contains information for 
%                                  one or more files. hFile is obtained by
%                                  a call to ns_OpenFile.  The result of 
%                                  ns_OpenFile, hFile, can be given as an 
%                                  argument in its entirety, or a subset of
%                                  the hFile array can be given.
%
% Return Values:
% nsFileInfo                       A cell array of ns_FILEINFO structures 
%                                  that receives the file information for 
%                                  each of the files referenced in the hFile 
%                                  handle. 
%
% struct ns_FILEINFO
% char FileType[32];               Human readable manufacturer�s 
%                                  file type descriptor.
%
% double EntityCount;              Number of entities in the data file.  
%                                  This number is used to enumerate all 
%                                  the entities in the data file from 1 
%                                  to EntityCount.
% 
% double TimeStampResolution       Minimum timestamp resolution in seconds.
% 
% double TimeSpan;                 Time span covered by the data file in 
%                                  seconds.
% 
% char AppName[64];                Information about the application 
%                                  that created the file.
% 
% double Time_Year;                Year.
% 
% double Time_Month;               Month (0-11; January = 0).
% 
% double Time_Day;                 Day of the month (1-31).
% 
% double Time_Hour;                Hour since midnight (0-23).
% 
% double Time_Min;                 Minute after the hour (0-59).
% 
% double Time_Sec;                 Seconds after the minute (0-59).
% 
% double Time_MilliSec;            Milliseconds after the second (0-1000).
% 
% char FileComment[256];           Comments embedded in the source file.
% 
% ns_RESULT          This function returns one of the following status 
%                    codes:
%
%   ns_OK              The file was successfully opened. 
%   ns_FILEERROR       File access or read error.
%   ns_BADFILE         Invalid file handle passed to function.
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
nsFileInfo = [];
if ~isstruct(hFile)
    nsFileInfo = [];
    ns_RESULT = 'ns_BADFILE';
    return
end
% initialize nsFileInfo structure
ns_RESULT = 'ns_OK';
nsFileInfo.FileType = sprintf('%s, ', hFile.FileInfo.Type);
nsFileInfo.FileType = nsFileInfo.FileType(1:end-2);
nsFileInfo.EntityCount = length(hFile.Entity);
nsFileInfo.AppName = [];
nsFileInfo.FileComment = [];
nsFileInfo.TimeSpan = hFile.TimeSpan/30000;
nsFileInfo.TimeStampResolution = min([hFile.FileInfo.Period])/30000;
nsFileInfo.Time_Year = [];
nsFileInfo.Time_Month = [];
nsFileInfo.Time_Day = [];
nsFileInfo.Time_Hour = [];
nsFileInfo.Time_Min = [];
nsFileInfo.Time_Sec = [];
nsFileInfo.Time_MilliSec = [];

for i = 1:length(hFile.FileInfo)
    out = fopen(hFile.FileInfo(i).FileID);
    if isempty(out)
        ns_RESULT = 'ns_FILEERROR';
        return
    end
    if strcmp(hFile.FileInfo(i).Type, 'nev')
        % seek to Time Origin header field
        fseek(hFile.FileInfo(i).FileID, 28, -1);
        Date = fread(hFile.FileInfo(i).FileID, 8, 'uint16');
        nsFileInfo.Time_Year = Date(1);
        nsFileInfo.Time_Month = Date(2);
        nsFileInfo.Time_Day = Date(4);
        nsFileInfo.Time_Hour = Date(5);
        nsFileInfo.Time_Min = Date(6);
        nsFileInfo.Time_Sec = Date(7);
        nsFileInfo.Time_MilliSec = Date(8);
        nsFileInfo.AppName =...
            deblank(fread(hFile.FileInfo(i).FileID, 32, '*char')');
        % In Ripple files (recorded with Trellis) the comment field is only
        % 252 characters and the NIP timestamp (uint32) for the recording
        % started is packed in the last four bytes.
        if ~isempty(strfind(nsFileInfo.AppName, 'Trellis'))
           nsFileInfo.FileComment =...
              deblank(fread(hFile.FileInfo(i).FileID, 252, '*char')');
           % With Ripple data get the NIP time
           nsFileInfo.NIPTime = fread(hFile.FileInfo(i).FileID, 1, 'uint32');
        else
           nsFileInfo.FileComment =...
              deblank(fread(hFile.FileInfo(i).FileID, 256, '*char')');
        end
    else % if not nev
        % check if nsx 2.2 file and that file comment is not yet written
        if isempty(nsFileInfo.FileComment)&&...
                (strcmp(hFile.FileInfo(i).FileTypeID, 'NEURALCD'))
            % seek to file comment
            fseek(hFile.FileInfo(i).FileID, 30, -1);
            % In Ripple files (recorded with Trellis) the comment field is only
            % 252 characters and the NIP timestamp (uint32) for the recording
            % started is packed in the last four bytes.
            %
            % We start first read the first 252 chars of the comment
            nsFileInfo.FileComment =...
               fread(hFile.FileInfo(i).FileID, 252, '*char')';

            if ~isempty(strfind(nsFileInfo.FileComment, 'Trellis'))
               % With Ripple data get the NIP time
               nsFileInfo.NIPTime = fread(hFile.FileInfo(i).FileID, 1, 'uint32');
               % NSX 2.2 doesn't have an explicit space for app name so
               % it's packed in the FileComment.  We'll repeat this in the
               % AppName
               nsFileInfo.FileComment = deblank(nsFileInfo.FileComment);
               nsFileInfo.AppName = nsFileInfo.FileComment;
            else
               nsFileInfo.FileComment = deblank([nsFileInfo.FileComment, ...
                  fread(hFile.FileInfo(i).FileID, 4, '*char')]);
            end
            % seek to Time Origin
            fseek(hFile.FileInfo(i).FileID, 8, 0);
            Date = fread(hFile.FileInfo(i).FileID, 8, 'uint16');
            nsFileInfo.Time_Year = Date(1);
            nsFileInfo.Time_Month = Date(2);
            nsFileInfo.Time_Day = Date(4);
            nsFileInfo.Time_Hour = Date(5);
            nsFileInfo.Time_Min = Date(6);
            nsFileInfo.Time_Sec = Date(7);
            nsFileInfo.Time_MilliSec = Date(8);
        end
    end
end
    