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

if isfield(res_, 'hClust')
    if ~isempty(res_.hClust.inconsistentFields()) && obj.hCfg.getOr('autoRecover', 0)
        flag = res_.hClust.recover(1); % recover inconsistent data if needed
        successAppend = 'You should look through your data and ensure everything is correct, then save it.';
        failureAppend = 'You will probably experience problems curating your data.';
        msg = '';
        switch flag
            case 2
                msg = sprintf('Non-contiguous spike table found and corrected. %s', successAppend);

            case 1
                msg = sprintf('Inconsistent fields found and corrected. %s', successAppend);

            case 0
                msg = sprintf('Clustering data in an inconsistent state and automatic recovery failed. Please post an issue on the GitHub issue tracker. %s', failureAppend);

            case -1
                msg = sprintf('Automatic recovery canceled by the user but the clustering data is still in an inconsistent state. %s', failureAppend);
        end

        if ~isempty(msg)
            jrclust.utils.qMsgBox(msg, 1, 1);
        end
    end

    % convert old-style history to new-style
    if isprop(obj.hCfg, 'histFile') && exist(obj.hCfg.histFile, 'file') == 2
        res_.hClust.history = obj.convertHistory();
        delete(obj.hCfg.histFile);
    end
end

obj.res = res_;
end

