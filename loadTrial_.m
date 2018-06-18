%--------------------------------------------------------------------------
% 7/26/17: Now accept csv file format
function vrTime_trial = loadTrial_(vcFile_trial)
    % import  trial time (in seconds)
    vrTime_trial = [];
    try
        if ~exist(vcFile_trial, 'file'), return; end

        if matchFileExt_(vcFile_trial, '.mat')
            Strial = load(vcFile_trial);
            csFields = fieldnames(Strial);
            vrTime_trial = Strial.(csFields{1});
            if isstruct(vrTime_trial)
                vrTime_trial = vrTime_trial.times;
            end
        elseif matchFileExt_(vcFile_trial, '.csv')
            vrTime_trial = csvread(vcFile_trial);
            if isrow(vrTime_trial), vrTime_trial = vrTime_trial(:); end
        end
    catch
        disperr_();
    end
end %func
