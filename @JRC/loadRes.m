function loadRes(obj)
%LOADRES Load results struct, storing it within obj.
%   Load data from hCfg.resFile, creating hClust if applicable.
if exist(obj.hCfg.resFile, 'file') ~= 2
    obj.error('%s does not exist', obj.hCfg.resFile);
end

try
    obj.hCfg.updateLog('loadRes', sprintf('Loading %s', obj.hCfg.resFile), 1, 0);
    res_ = load(obj.hCfg.resFile);
    obj.hCfg.updateLog('loadRes', sprintf('Finished loading %s', obj.hCfg.resFile), 0, 1);
catch ME
    warning('Failed to load %s: %s', ME.message);
    return;
end

%% check for Clustering fields, create hClust if found
if isfield(res_, 'spikeTemplates') % create a new TemplateClustering
    hClust = jrclust.sort.TemplateClustering(obj.hCfg);
    fieldNames = fieldnames(res_);

    md = ?jrclust.sort.TemplateClustering;
    pl = md.PropertyList;
    for i = 1:numel(fieldNames)
        fn = fieldNames{i};
        propMetadata = pl(strcmp(fn, {pl.Name}));
        if isprop(hClust, fn) && ~((propMetadata.Dependent && isempty(propMetadata.SetMethod)))
            hClust.(fn) = res_.(fn);
        end
    end

    res_.hClust = hClust;
elseif isfield(res_, 'spikeClusters')
    hClust = jrclust.sort.DensityPeakClustering(obj.hCfg);
    fieldNames = fieldnames(res_);

    md = ?jrclust.sort.DensityPeakClustering;
    pl = md.PropertyList;
    for i = 1:numel(fieldNames)
        fn = fieldNames{i};
        propMetadata = pl(strcmp(fn, {pl.Name}));
        if isprop(hClust, fn) && ~((propMetadata.Dependent && isempty(propMetadata.SetMethod)))
            try
                hClust.(fn) = res_.(fn);
            catch
            end
        end
    end

    res_.hClust = hClust;
end

obj.res = res_;
end

