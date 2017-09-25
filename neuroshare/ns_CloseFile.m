function ns_RESULT = ns_CloseFile(hFile)
% ns_CloseFile - Closes a neural data file
% 
% Usage:
% ns_RESULT = ns_CloseFile(hFile)
% 
% Description:
% Closes previously opened file(s) specified by the file handle hFile.
%
% Parameters:
% 
% hFile                            A handle that contains information for 
%                                  one or more files. hFile is obtained by
%                                  a call to ns_OpenFile.  The result of 
%                                  ns_OpenFile, hFile, can be given as an 
%                                  argument in its entirety, or a subset of
%                                  the hFile array can be given.
%
% Return Values:
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

% Ensure hFile looks like file handle returned from ns_OpenFile
% sometimes, hFile get set to a scalar or [] and this function fails.
if ~isstruct(hFile)
    ns_RESULT = 'ns_BADFILE';
    return
end

ns_RESULT = 'ns_OK';

fids     = [hFile.FileInfo.FileID];
errCount = 0;
for i=1:length(fids)
    errCount = errCount + fclose(fids(i));
end

if errCount
    ns_RESULT = 'ns_BADFILE';
end

% NOTE: [AMW-2012/04/12] 
% For some reason this does not to work - it leaves file handles open
%
% try
%     arrayfun(@fclose, [hFile.FileInfo.FileID]);
% catch
%     ns_RESULT = 'ns_BADFILE';
end

