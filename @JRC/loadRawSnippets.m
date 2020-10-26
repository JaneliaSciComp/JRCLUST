function spikesRaw = loadRawSnippets(obj)
%LOADRAWSNIPPETS Load raw spike snippets from disk.
spikesRaw = [];

if isfield(obj.res, 'rawShape') && ~isempty(obj.res.rawShape) && exist(obj.hCfg.rawFile, 'file') == 2
    try
        obj.hCfg.updateLog('loadRaw', sprintf('Loading %s', obj.hCfg.rawFile), 1, 0);
        spikesRaw = readBin(obj.hCfg.rawFile, obj.res.rawShape, '*int16');
        obj.hCfg.updateLog('loadRaw', sprintf('Finished loading %s', obj.hCfg.rawFile), 0, 1);
    catch ME
        obj.hCfg.updateLog('loadRaw', sprintf('Failed to load %s: %s', obj.hCfg.rawFile, ME.message), 0, 1);
    end
end
end %fun

