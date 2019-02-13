function saveRes(obj, forceOverwrite)
    %SAVERES Save results struct
    if isempty(obj.res)
        return;
    end
    if nargin < 2
        forceOverwrite = 0;
    end
    forceOverwrite = forceOverwrite | obj.hCfg.getOr('testRun', 0);

    sessionName = obj.hCfg.sessionName;
    filename = fullfile(obj.hCfg.outputDir, [sessionName '_res.mat']);
     % don't overwrite unless explicitly permitted (or testing)
    if exist(filename, 'file') && ~forceOverwrite
        question = sprintf('%s already exists. Overwrite?', filename);
        dlgAns = questdlg(question, 'Confirm overwrite', 'No');
        if isempty(dlgAns) || ismember(dlgAns, {'No', 'Cancel'})
            return;
        end
    end

    % save everything else (don't save spikesRaw, spikesFilt,
    % spikeFeatures inside hClust)
    if isfield(obj.res, 'hClust') && ~isempty(obj.res.hClust.spikesRaw)
        doRestore = 1;

        spikesRaw = obj.res.hClust.spikesRaw;
        obj.res.hClust.spikesRaw = [];

        spikesFilt = obj.res.hClust.spikesFilt;
        obj.res.hClust.spikesFilt = [];

        spikeFeatures = obj.res.hClust.spikeFeatures;
        obj.res.hClust.spikeFeatures = [];
    else
        doRestore = 0;
    end

    jrclust.utils.saveStruct(obj.res, filename);

    if doRestore % restore spikesRaw, spikesFilt, spikeFeatures to hClust
        obj.res.hClust.spikesRaw = spikesRaw;
        obj.res.hClust.spikesFilt = spikesFilt;
        obj.res.hClust.spikeFeatures = spikeFeatures;
    end

    if obj.hCfg.verbose
        fprintf('Saved results to %s\n', filename);
    end
end