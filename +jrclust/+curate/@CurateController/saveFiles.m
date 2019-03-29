function success = saveFiles(obj)
    %SAVEFILES Save files to disk
    dlgAns = questdlg(['Save clustering to ', obj.hCfg.resFile, ' ?'], 'Confirmation', 'Yes');
    if strcmp(dlgAns, 'Yes')
        obj.cRes.curatedOn = now();

        res_ = jrclust.utils.mergeStructs(obj.res, obj.cRes);
        res_ = rmfield(res_, 'hClust');

        hClust = obj.hClust;

        res_.initialClustering = hClust.initialClustering;
        res_.spikeClusters = hClust.spikeClusters;

        % fieldnames contained in dRes or sRes
        dsFields = union(fieldnames(hClust.dRes), fieldnames(hClust.sRes));
        % fieldnames from hClust which are not in dRes or sRes
        hClustOnly = setdiff(fieldnames(hClust), dsFields);
        for i = 1:numel(hClustOnly)
            fn = hClustOnly{i};
            % hClust fields take precedence
            res_.(fn) = hClust.(fn);
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

        hBox = jrclust.utils.qMsgBox('Saving... (this closes automatically)', 0, 1);
        jrclust.utils.saveStruct(res_, obj.hCfg.resFile);

        success = 1;
        jrclust.utils.tryClose(hBox);
    elseif strcmp(dlgAns, 'No')
        success = 1;
    else
        success = 0;
    end
end