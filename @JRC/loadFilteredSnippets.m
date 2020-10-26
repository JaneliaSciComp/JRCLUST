function spikesFilt = loadFilteredSnippets(obj)
%LOADFILTEREDSNIPPETS Load filtered spike snippets from disk.
spikesFilt = [];

if isfield(obj.res, 'filtShape') && ~isempty(obj.res.filtShape) && exist(obj.hCfg.filtFile, 'file') == 2
    try
        obj.hCfg.updateLog('loadFilt', sprintf('Loading %s', obj.hCfg.filtFile), 1, 0);
        spikesFilt = readBin(obj.hCfg.filtFile, obj.res.filtShape, '*int16');
        obj.hCfg.updateLog('loadFilt', sprintf('Finished loading %s', obj.hCfg.filtFile), 0, 1);
    catch ME
        obj.hCfg.updateLog('loadFilt', sprintf('Failed to load %s: %s', obj.hCfg.filtFile, ME.message), 0, 1);
    end
end
end %fun

