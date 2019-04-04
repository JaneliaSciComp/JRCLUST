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

    if isfield(obj.res, 'hClust')
        hClust = obj.res.hClust;

        res_ = rmfield(obj.res, 'hClust');
        res_.initialClustering = hClust.initialClustering;
        res_.spikeClusters = hClust.spikeClusters;

        % fieldnames contained in dRes or sRes
        fieldNames = fieldnames(hClust);
        for i = 1:numel(fieldNames)
            fn = fieldNames{i};
            res_.(fn) = hClust.(fn);
        end
    else
        res_ = obj.res;
    end

    % don't save these fields
    if isfield(res_, 'spikesRaw')
        res_ = rmfield(res_, 'spikesRaw');
    end
    if isfield(res_, 'spikesRawVolt')
        res_ = rmfield(res_, 'spikesRawVolt');
    end
    if isfield(res_, 'spikesFilt')
        res_ = rmfield(res_, 'spikesFilt');
    end
    if isfield(res_, 'spikesFiltVolt')
        res_ = rmfield(res_, 'spikesFiltVolt');
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
    if isfield(res_, 'nClusters')
        res_ = rmfield(res_, 'nClusters');
    end
    if isfield(res_, 'nEdits')
        res_ = rmfield(res_, 'nEdits');
    end
    if isfield(res_, 'unitFields')
        res_ = rmfield(res_, 'unitFields');
    end
    if isfield(res_, 'hRecs')
        res_ = rmfield(res_, 'hRecs');
    end

    jrclust.utils.saveStruct(res_, obj.hCfg.resFile);

    obj.hCfg.updateLog('saveRes', sprintf('Results saved to %s', obj.hCfg.resFile), 0, 1);
end
