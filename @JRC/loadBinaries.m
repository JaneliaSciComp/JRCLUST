function loadBinaries(obj)
%LOADBINARIES Load binary files from disk.
%   Store raw and filtered spike snippets and clustering features in res
%   and hClust (if it exists).

% load raw spike snippets
obj.res.spikesRaw = obj.loadRawSnippets();

% load filtered spike snippets
obj.res.spikesFilt = obj.loadFilteredSnippets();
    
% load features
obj.res.spikeFeatures = obj.loadSpikeFeatures();

if isfield(obj.res, 'hClust')
    obj.hClust.spikesRaw = obj.res.spikesRaw;
    obj.hClust.spikesFilt = obj.res.spikesFilt;
    obj.hClust.spikeFeatures = obj.res.spikeFeatures;
end
end %fun