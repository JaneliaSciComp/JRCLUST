function autoMerge(obj)
    %AUTOMERGE Automatically merge units based on their similarity scores
    % snr_thresh = jrclust.utils.inputdlgNum('SNR threshold: ', 'Auto-deletion based on SNR', 10); % also ask about # spikes/unit (or firing rate) @TODO
    if obj.isWorking
        return;
    end
    obj.isWorking = 1;

    dlgAns = inputdlg('Waveform correlation threshold (0-1):', 'Auto-merge based on waveform threshold', 1, {num2str(obj.hCfg.maxUnitSim)});

    % parse user input
    if isempty(dlgAns)
        return;
    end

    mwc = str2double(dlgAns{1});
    if isnan(mwc) || mwc <= 0 || mwc > 1
        jrclust.utils.qMsgBox('Invalid criteria.');
        return;
    end

    % auto merge
    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(1);
    end

    nClustersOld = obj.hClust.nClusters;

    hBox = jrclust.utils.qMsgBox('Merging...', 0, 1);
    if obj.hClust.autoMerge(mwc) % success; replot
        jrclust.utils.tryClose(hBox);
        obj.updateFigWav();
        obj.updateFigRD(); % centers changed, need replotting
        obj.updateFigSim();
        obj.updateSelect(obj.selected);

        jrclust.utils.qMsgBox(sprintf('Merged %d clusters >%0.2f maxUnitSim.', nClustersOld - obj.hClust.nClusters, mwc));
    else
        jrclust.utils.tryClose(hBox);
        jrclust.utils.qMsgBox('Auto merge failed.');
    end

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(0);
    end

    obj.isWorking = 0;
end