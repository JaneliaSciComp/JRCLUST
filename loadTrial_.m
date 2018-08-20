%--------------------------------------------------------------------------
% 7/26/17: Now accept csv file format
function vrTime_trial = loadTrial_(trialFile)
    % import  trial time (in seconds)
    vrTime_trial = [];
    try
        if ~exist(trialFile, 'file'), return; end

        if matchFileExt_(trialFile, '.mat')
            Strial = load(trialFile);
            csFields = fieldnames(Strial);
            vrTime_trial = Strial.(csFields{1});
            if isstruct(vrTime_trial)
                vrTime_trial = vrTime_trial.times;
            end
        elseif matchFileExt_(trialFile, '.csv')
            vrTime_trial = csvread(trialFile);
            if isrow(vrTime_trial), vrTime_trial = vrTime_trial(:); end
        end
    catch
        disperr_();
    end
end %func
