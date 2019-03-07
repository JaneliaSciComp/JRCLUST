function autoMerge(obj, maxUnitSim)
    %AUTOMERGE Automatically merge units based on their similarity scores
    % snr_thresh = jrclust.utils.inputdlgNum('SNR threshold: ', 'Auto-deletion based on SNR', 10); % also ask about # spikes/unit (or firing rate) @TODO
    if obj.isWorking
        jrclust.utils.qMsgBox('An operation is in progress.');
        return;
    elseif ~isa(obj.hClust, 'jrclust.sort.DensityPeakClustering')
        jrclust.utils.qMsgBox('Operation not supported for this type of clustering');
        return;
    end

    if nargin < 2
        dlgAns = inputdlg('Waveform correlation threshold (0-1):', 'Auto-merge based on waveform threshold', 1, {num2str(obj.hCfg.maxUnitSim)});
        % parse user input
        if isempty(dlgAns)
            return;
        end

        maxUnitSim = str2double(dlgAns{1});
        if isnan(maxUnitSim) || maxUnitSim <= 0 || maxUnitSim > 1
            jrclust.utils.qMsgBox('Invalid criteria.');
            return;
        end
    end

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(1);
    end

    nClustersOld = obj.hClust.nClusters;
    hBox = jrclust.utils.qMsgBox('Merging... (this closes automatically)', 0, 1);

    obj.isWorking = 1;
    try
        obj.hClust.hCfg.setTemporaryParams('maxUnitSim', maxUnitSim);
        if obj.hClust.autoMerge() % success; replot
            jrclust.utils.tryClose(hBox);

            obj.isWorking = 0; % in case updateSelect needs to zoom
            
            obj.updateFigWav();
            obj.updateFigRD(); % centers changed, need replotting
            obj.updateFigSim();
            obj.updateSelect(obj.selected);

            jrclust.utils.qMsgBox(sprintf('Merged %d clusters >%0.2f maxUnitSim.', nClustersOld - obj.hClust.nClusters, maxUnitSim));
        else
            jrclust.utils.tryClose(hBox);
            jrclust.utils.qMsgBox('Auto merge failed.');
        end

        obj.hClust.hCfg.resetTemporaryParams('maxUnitSim');

        if obj.hasFig('FigWav')
            hFigWav = obj.hFigs('FigWav');
            hFigWav.wait(0);
        end
    catch ME
        jrclust.utils.tryClose(hBox);
        rethrow(ME)
        warning('Failed to merge: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed.');
    end

    obj.isWorking = 0;
end