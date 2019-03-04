function recomputeWaveforms(obj, computeAll)
    %RECOMPUTEWAVEFORMS Recompute mean waveforms for a cluster or all
    %clusters

    hBox = jrclust.utils.qMsgBox('Recomputing... (this closes automatically)', 0, 1);

    obj.isWorking = 1;
    try
        if nargin < 2 || computeAll == 0
            obj.hClust.updateWaveforms(obj.selected);
        else
            obj.hClust.updateWaveforms([]);
        end
    catch ME
        warning('Operation failed: %s', ME.message);
        jrclust.utils.qMsgBox('Operation failed');
    end

    jrclust.utils.tryClose(hBox);
    obj.isWorking = 0;

    obj.updateFigWav();
    obj.updateFigSim();
    obj.updateSelect(obj.selected);
end

