function plot_ripple_analog(hFile, port, headstage, channel)
%PLOT_RIPPLE_ANALOG Plot analog data from Grapevine port ids
% This function is a wrapper for PLOT_ANALOG, but allows for parameters
% that translate easily to Ripple Grapevine NIP ports and headstages.
% 
% Usage: ns_RESULT = plot_ripple_analog(hFile, port, headstage, channel)
%   Parameters:
%       hFile - Neuroshare file handle (returned from ns_OpenFile)
%       port - single character that should be 'a' - 'd'
%       headstage - index of headstage used
%       channel - index of channel wanted relative to the headstage.
% See also PLOT_ANALOG

[ns_RESULT, entityIDs] = ...
    get_ripple_entity(hFile, port, headstage, channel);
if ~strcmp(ns_RESULT, 'ns_OK')
    fprintf(2, '%s error\n', ns_RESULT);
    return
end 
analog_entities = entityIDs(strcmp({entityIDs(:).type}, 'Analog'));
analog_length = length(analog_entities);
if analog_length > 2
    fprintf(2, 'error: too many analog entities\n');
    return
end
if analog_length == 0
    fprintf(2, 'error: cannot find analog entities\n');
    return
end

file_indexes = [hFile.Entity([analog_entities(:).index]).FileType];
periods = [hFile.FileInfo(file_indexes).Period];
%[analog_entities(:).index];
if periods(1) == 30
    entity = analog_entities(1).index;
elseif length(periods) == 2 && periods(2) == 1
    entity = analog_entities(2).index;
else
    fprintf(2, 'error: cannot find 1 kS data\n');
    return
end    
plot_analog(hFile, entity);

