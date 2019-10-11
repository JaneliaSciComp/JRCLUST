function trialTimes = loadTrialFile(trialFile)
    %LOADTRIALFILE Import trial times (in seconds)
    trialTimes = [];

    try
        [~, ~, ext] = fileparts(trialFile);

        if strcmpi(ext, '.mat')
            trialData = load(trialFile);
            fieldNames = fieldnames(trialData);

            trialTimes = trialData.(fieldNames{1});
            if isstruct(trialTimes)
                trialTimes = trialTimes.times;
            end
        elseif strcmpi(ext, '.csv')
            trialTimes = csvread(trialFile);
            if isrow(trialTimes)
                trialTimes = trialTimes(:);
            end
        end
    catch ME
        warning('Could not load trialFile %s: %s', trialFile, ME.message);
    end
end