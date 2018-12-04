function [dRes, sRes, cRes] = doLoadFiles(obj)
    %DOLOADFILES Load results structs
    [dRes, sRes, cRes] = deal([]);

    [~, sessionName, ~] = fileparts(obj.hCfg.configFile);
    if nargout >= 1
        filename = fullfile(obj.hCfg.outputDir, [sessionName '_detect.mat']);
        try
            dRes = load(filename);
        catch ME
            warning(ME.identifier, 'detect data not loaded: %s', ME.message);
        end
    end
    if nargout >= 2
        filename = fullfile(obj.hCfg.outputDir, [sessionName '_sort.mat']);
        try
            sRes = load(filename);
        catch ME
            warning(ME.identifier, 'sort data not loaded: %s', ME.message);
        end
    end
    if nargout >= 3
        filename = fullfile(obj.hCfg.outputDir, [sessionName '_curate.mat']);
        try
            cRes = load(filename);
        catch ME
            warning(ME.identifier, 'curate data not loaded: %s', ME.message);
        end
    end
end

