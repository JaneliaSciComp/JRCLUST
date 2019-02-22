function saveRes(obj, forceOverwrite)
    %SAVERES Save results struct
    if isempty(obj.res)
        return;
    end
    if nargin < 2
        forceOverwrite = 0;
    end
    forceOverwrite = forceOverwrite | obj.hCfg.getOr('testRun', 0);

    % don't overwrite unless explicitly permitted (or testing)
    if exist(obj.hCfg.resFile, 'file') && ~forceOverwrite
        question = sprintf('%s already exists. Overwrite?', obj.hCfg.resFile);
        dlgAns = questdlg(question, 'Confirm overwrite', 'No');
        if isempty(dlgAns) || ismember(dlgAns, {'No', 'Cancel'})
            return;
        end
    end

    obj.hCfg.updateLog('saveRes', sprintf('Saving results to %s', obj.hCfg.resFile), 1, 0);
    % save everything else (don't save spikesRaw, spikesFilt,
    % spikeFeatures inside hClust)
    restoreFields = struct(); % restore these fields to hClust after saving
    if isfield(obj.res, 'hClust')
        restoreFields.dRes = obj.res.hClust.dRes;
        obj.res.hClust.dRes = [];

        restoreFields.sRes = obj.res.hClust.sRes;
        obj.res.hClust.sRes = [];

        % fieldnames contained in dRes and sRes
        dsFields = union(fieldnames(restoreFields.dRes), fieldnames(restoreFields.sRes));
        % fieldnames from hClust which are not in dRes or sRes
        hClustOnly = setdiff(fieldnames(obj.res.hClust), dsFields);
        % any redundancies with res?
        intersectFields = intersect(fieldnames(obj.res), hClustOnly);
        for i = 1:numel(intersectFields)
            fn = intersectFields{i};
            if jrclust.utils.isEqual(obj.res.(fn), obj.res.hClust.(fn))
                restoreFields.(fn) = obj.res.hClust.(fn);
                obj.res.hClust.(fn) = [];
            end
        end
    end

    % don't save these large fields
    res_ = obj.res;
    if isfield(res_, 'spikesRaw')
        res_ = rmfield(res_, 'spikesRaw');
    end
    if isfield(res_, 'spikesFilt')
        res_ = rmfield(res_, 'spikesFilt');
    end
    if isfield(res_, 'spikesFilt2')
        res_ = rmfield(res_, 'spikesFilt2');
    end
    if isfield(res_, 'spikesFilt3')
        res_ = rmfield(res_, 'spikesFilt3');
    end
    if isfield(res_, 'spikeFeatures')
        res_ = rmfield(res_, 'spikeFeatures');
    end
    if isfield(res_, 'hRecs')
        res_ = rmfield(res_, 'hRecs');
    end

    jrclust.utils.saveStruct(res_, obj.hCfg.resFile);

    restoreFieldNames = fieldnames(restoreFields);
    for i = 1:numel(restoreFieldNames)
        fn = restoreFieldNames{i};
        obj.res.hClust.(fn) = restoreFields.(fn);
    end

    obj.hCfg.updateLog('saveRes', sprintf('Results saved to %s', obj.hCfg.resFile), 0, 1);
end