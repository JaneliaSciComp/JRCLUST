function [needsConversion, needsDeletion] = checkHistory(obj)
%CHECKHISTORY Check for old-style history file and determine if history
%needs conversion.
%   History in hClust 
needsConversion = false;
needsDeletion = false;

if exist(obj.hCfg.histFile, 'file') == 2
    needsDeletion = true;
end

% nothing to convert
if ~needsDeletion
    return;
end

% nothing to convert if file is empty
d = dir(obj.hCfg.histFile);
if d.bytes == 0
    return;
end

% no need to convert a history file when there's no Clustering
if ~isfield(obj.res, 'hClust')
    return;
end

% no need to convert a history file if the history field is already in the
% correct format
history = obj.hClust.history;
if isstruct(history) && ...
        isfield(history, 'optype') && iscell(history.optype) && ... % optype field exists and is a cell
        isfield(history, 'message') && iscell(history.message) && ... % message field exists and is a cell
        isfield(history, 'indices') && iscell(history.indices) && ... % indices field exists and is a cell
        numel(history.optype) == numel(history.message) && ... % optype, message, and indices fields all have the same length
        numel(history.message) == numel(history.indices)
    return;
end

% still here? we need to convert that history
dlgAns = questdlg('Your history needs to be converted from the old style to the new. Please confirm.', ...
    'Confirm Convert History', 'OK', 'Cancel', 'OK');

if strcmp(dlgAns, 'OK')
    needsConversion = true;
end

end %fun

