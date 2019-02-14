function success = saveFiles(obj)
    %SAVEFILES Save files to disk
    dlgAns = questdlg(['Save clustering to ', obj.hCfg.resFile, ' ?'], 'Confirmation', 'Yes');
    if strcmp(dlgAns, 'Yes')
        obj.cRes.curatedOn = now();

        % don't save these to _res.mat
        spikesRaw = obj.hClust.spikesRaw;
        obj.hClust.spikesRaw = [];

        spikesFilt = obj.hClust.spikesFilt;
        obj.hClust.spikesFilt = [];

        spikeFeatures = obj.hClust.spikeFeatures;
        obj.hClust.spikeFeatures = [];

        hMsg = jrclust.utils.qMsgBox('Saving... (this closes automatically)');

        % save hClust and curatedOn
        cRes = obj.cRes; %#ok<NASGU>
        save(obj.hCfg.resFile, '-struct', 'cRes', '-append');

        % restore these values to hClust
        obj.hClust.spikesRaw = spikesRaw;
        obj.hClust.spikesFilt = spikesFilt;
        obj.hClust.spikeFeatures = spikeFeatures;

        success = 1;
        jrclust.utils.tryClose(hMsg);
    elseif strcmp(dlgAns, 'No')
        success = 1;
    else
        success = 0;
    end
end