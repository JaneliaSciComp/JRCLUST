function success = saveFiles(obj)
    %SAVEFILES Save files to disk
    dlgAns = questdlg(['Save clustering to ', obj.hCfg.resFile, ' ?'], 'Confirmation', 'Yes');
    if strcmp(dlgAns, 'Yes')
        obj.cRes.curatedOn = now();

        res = jrclust.utils.mergeStructs(obj.res, obj.cRes);

        % don't save this stuff twice
        restoreFields = struct(); % restore these fields to hClust after saving
        restoreFields.dRes = res.hClust.dRes;
        res.hClust.dRes = [];

        restoreFields.sRes = res.hClust.sRes;
        res.hClust.sRes = [];

        hMsg = jrclust.utils.qMsgBox('Saving... (this closes automatically)', 0, 1);
        jrclust.utils.saveStruct(res, obj.hCfg.resFile);

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