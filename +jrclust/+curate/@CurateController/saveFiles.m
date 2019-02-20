function success = saveFiles(obj)
    %SAVEFILES Save files to disk
    dlgAns = questdlg(['Save clustering to ', obj.hCfg.resFile, ' ?'], 'Confirmation', 'Yes');
    if strcmp(dlgAns, 'Yes')
        obj.cRes.curatedOn = now();

        res_ = jrclust.utils.mergeStructs(obj.res, obj.cRes);
        % don't save these large fields
        if isfield(res_, 'spikesRaw')
            res_ = rmfield(res_, 'spikesRaw');
        end
        if isfield(res_, 'spikesFilt')
            res_ = rmfield(res_, 'spikesFilt');
        end
        if isfield(res_, 'spikeFeatures')
            res_ = rmfield(res_, 'spikeFeatures');
        end

        % don't save this stuff twice
        restoreFields = struct(); % restore these fields to hClust after saving
        restoreFields.dRes = res_.hClust.dRes;
        res_.hClust.dRes = [];

        restoreFields.sRes = res_.hClust.sRes;
        res_.hClust.sRes = [];

        hMsg = jrclust.utils.qMsgBox('Saving... (this closes automatically)', 0, 1);
        jrclust.utils.saveStruct(res_, obj.hCfg.resFile);

        % restore fields to carry on
        if ~obj.isEnding
            restoreFieldNames = fieldnames(restoreFields);
            for i = 1:numel(restoreFieldNames)
                fn = restoreFieldNames{i};
                obj.res.hClust.(fn) = restoreFields.(fn);
            end
        end

        success = 1;
        jrclust.utils.tryClose(hMsg);
    elseif strcmp(dlgAns, 'No')
        success = 1;
    else
        success = 0;
    end
end