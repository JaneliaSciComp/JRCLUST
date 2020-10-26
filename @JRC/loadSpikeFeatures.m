function spikeFeatures = loadSpikeFeatures(obj)
%LOADSPIKEFEATURES Load clustering features from disk.
spikeFeatures = [];

if isfield(obj.res, 'featuresShape') && ~isempty(obj.res.featuresShape) && exist(obj.hCfg.featuresFile, 'file') == 2
    try
        obj.hCfg.updateLog('loadFeatures', sprintf('Loading %s', obj.hCfg.featuresFile), 1, 0);
        spikeFeatures = readBin(obj.hCfg.featuresFile, obj.res.featuresShape, '*single');
        obj.hCfg.updateLog('loadFeatures', sprintf('Finished loading %s', obj.hCfg.featuresFile), 0, 1);
    catch ME
        obj.hCfg.updateLog('loadFeatures', sprintf('Failed to load %s: %s', obj.hCfg.featuresFile, ME.message), 0, 1);
    end
end
end %fun

