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
            hClust.(fn) = res_.(fn);
        end
    end

    res_.hClust = hClust;
end

%% convert history if necessary
if isfield(res_, 'hClust')
    % convert old-style history to new-style
    if isprop(obj.hCfg, 'histFile') && exist(obj.hCfg.histFile, 'file') == 2
        safeToDelete = 0;
        try
            history = obj.convertHistory();
            safeToDelete = 1;
        catch ME
            warning('Failed to convert history: %s', ME.message);
        end

        % delete superfluous history file
        if safeToDelete
            try
                delete(obj.hCfg.histFile);
            catch ME
                warning('Failed to delete %s: %s', ME.message);
            end
        end

        res_.hClust.history = history;
        delete(obj.hCfg.histFile);
    end
end

obj.res = res_;
end

