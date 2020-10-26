function loadBinaries(obj)
%LOADBINARIES Load binary files from disk.
%   Store raw and filtered spike snippets and clustering features.

if isfield(obj.res, 'hClust')
    % load raw spike snippets
    obj.hClust.spikesRaw = obj.loadRawSnippets();

    % load filtered spike snippets
    obj.hClust.spikesFilt = obj.loadFilteredSnippets();
    
    % load features
    obj.hClust.spikeFeatures = obj.loadSpikeFeatures();
end
end