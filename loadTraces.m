%--------------------------------------------------------------------------
function S0 = loadTraces(P, spikeTimes0, spikeSites0)
    S0 = [];

    spikeTimes0 = spikeTimes0(:);
    spikeSites0 = spikeSites0(:);

    % regularize list of files to load
    if ~isfield(P, 'multiFilenames') || isempty(P.multiFilenames)
        if ~fileExists(P.vcFile)
            P.vcFile = replaceDir(P.vcFile, P.paramFile);
        end

        filenames = {P.vcFile};
    else
        filenames = filter_files_(P.multiFilenames);
        if isempty(filenames)
            P.multiFilenames = replaceDir(P.multiFilenames, P.paramFile);
            filenames = filter_files_(P.multiFilenames);
        end
    end

    % load site-wise thresholds, if applicable
    if ~isempty(get_(P, 'thresholdFile'))
        try
            S_thresh = load(P.thresholdFile);
            siteThresholds = S_thresh.siteThresholds;
            setUserData(siteThresholds);
            fprintf('Loaded %s\n', P.thresholdFile);
        catch
            disperr_('thresholdFile load error');
        end
    end

    if isempty(filenames)
        error('No binary files found.');
    end

    nFiles = numel(filenames);
end % fund
