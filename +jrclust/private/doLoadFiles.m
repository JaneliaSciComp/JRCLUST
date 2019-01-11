function res = doLoadFiles(hCfg)
    %DOLOADFILES Load results struct
    res = struct();

    [~, sessionName, ~] = fileparts(hCfg.configFile);
    filename = fullfile(hCfg.outputDir, [sessionName '_res.mat']);
    if ~isfile(filename)
        warning('%s does not exist', filename);
        return;
    end

    try
        res = load(filename);
    catch ME
        warning(ME.identifier, 'failed to load %s: %s', ME.message);
    end
end

