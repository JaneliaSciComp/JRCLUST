function [ns_RESULT, LastError] = ns_GetLastErrorMsg
% ns_GetLastErrorMsg - Retrieves the extended error message for the last 
%                      error
% Usage:
% 
% [ns_RESULT, LastError] = ns_GetLastErrorMsg
% Description:
% 
% Returns the last error message in text form to LastError. This function 
% should be called immediately following a function whose return value 
% indicates that an error occurred. The maximum size of the error message 
% text is 256 characters.
%
% Return Values:
% 
% LastError         Variable to receive the extended last error message.
% 
% ns_RESULT         This function returns ns_OK if the file is successfully
%                   opened.

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
ns_RESULT = 'ns_OK';
LastError = lasterror.message;
end