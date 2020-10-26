function fieldNames = getSaveFields(obj)
%GETSAVEFIELDS Return a list of fields to save to or load from res file.
%   Relevant fields are from res struct and, where applicable, hClust.
if isempty(obj.res)
    fieldNames = {};
    return;
end

fieldNames = setdiff(fieldnames(obj.res), ...
    {'spikesRaw', 'spikesRawVolt', 'spikesFilt', 'spikesFiltVolt', ...
     'spikesFilt2', 'spikesFilt3', 'spikeFeatures', 'hClust', 'hRecs'});
end

