function doSaveFiles(hJRC)
    %DOSAVEFILES Save results structs and binary files to disk
    [~, sessionName, ~] = fileparts(hJRC.hCfg.configFile);

    if ~isempty(hJRC.dRes)
        filename = fullfile(hJRC.hCfg.outputDir, [sessionName '_detect.mat']);
        if hJRC.hCfg.verbose
            fprintf('Saving detection results to %s\n', filename);
        end

        jrclust.utils.saveStruct(hJRC.dRes, filename);
    end

    if ~isempty(hJRC.sRes)
        filename = fullfile(hJRC.hCfg.outputDir, [sessionName '_sort.mat']);
        if hJRC.hCfg.verbose
            fprintf('Saving sorting results to %s\n', filename);
        end

        sRes = hJRC.sRes;
        sRes.hClust = hJRC.hClust;
        sRes.hCfg = hJRC.hCfg;
        jrclust.utils.saveStruct(sRes, filename);
    end
end

