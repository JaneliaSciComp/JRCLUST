function trialTimes = loadTrialFile(trialFile)
    %LOADTRIALFILE Import trial times (in seconds)
    trialTimes = [];
    if ~exist(trialFile,'file')
        warning('%s does not exist.',trialFile);
        return
    end
    try
        [~, ~, ext] = fileparts(trialFile);
        if strcmpi(ext, '.mat') % use the first variable if there are multiple variables in the file. If that variable is a struct, use the field "times"
            trialData = load(trialFile);
            fieldNames = fieldnames(trialData);
            trialTimes = trialData.(fieldNames{1});
            if isstruct(trialTimes)
                if isfield(trialTimes,'times')
                    trialTimes = trialTimes.times;                
                else
                    warning('No field "times" in structure variable %s found in: %s',fieldNames{1},trialFile);
                    return
                end
            end
        elseif strcmpi(ext, '.csv')
            trialTimes = csvread(trialFile);
            if isrow(trialTimes)
                trialTimes = trialTimes(:);
            end
        end
        if isempty(trialTimes) || ~isvector(trialTimes)
            warning('trialTimes found is either empty or not a vector!');
            return
        end
    catch ME
        warning('Could not load trial information %s: %s', trialFile, ME.message);
    end
end