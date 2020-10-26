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

    % save fields from obj (JRC)
    res_ = struct();
    fieldNames = obj.getSaveFields();
    for i = 1:numel(fieldNames)
        fn = fieldNames{i};
        res_.(fn) = obj.res.(fn);
    end

    % save fields from obj.res.hClust, if applicable
    if isfield(obj.res, 'hClust')
        try
            fieldNames = obj.hClust.getSaveFields(); % obj.hClust -> obj.res.hClust
            for i = 1:numel(fieldNames)
                fn = fieldNames{i};
                res_.(fn) = obj.hClust.(fn);
            end
        catch ME
            warning(ME.message);
        end
    end

    jrclust.utils.saveStruct(res_, obj.hCfg.resFile);

    obj.hCfg.updateLog('saveRes', sprintf('Results saved to %s', obj.hCfg.resFile), 0, 1);
end
