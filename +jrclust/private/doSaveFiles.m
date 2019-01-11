function doSaveFiles(res, hCfg)
    %DOSAVEFILES Save results structs and binary files to disk
    [~, sessionName, ~] = fileparts(hCfg.configFile);

    if isempty(res)
        return;
    end

    filename = fullfile(hCfg.outputDir, [sessionName '_res.mat']);
    if isfile(filename) % don't overwrite without explicit permission
        question = sprintf('%s already exists. Overwrite?', filename);
        dlgAns = questdlg(question, 'Confirm overwrite', 'No');
        if isempty(dlgAns) || ismember(dlgAns, {'No', 'Cancel'})
            return;
        end
    end

    if hCfg.verbose
        fprintf('Saving results to %s\n', filename);
    end
    jrclust.utils.saveStruct(res, filename);
end

