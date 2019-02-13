function toggleRaw(obj, hMenu)
    %TOGGLERAW Toggle raw waveform display
    showRaw = ~obj.hCfg.showRaw;

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(1);
    end

    if showRaw && isempty(obj.hClust.meanWfGlobalRaw)
        obj.hClust.computeMeanWaveforms([], 1);
    end

    set(hMenu, 'Checked', jrclust.utils.ifEq(showRaw, 'on', 'off'));
    obj.hCfg.showRaw = showRaw;

    % replot
    obj.updateFigWav();
    obj.updateFigSim();
    obj.updateSelect(obj.selected);

    if obj.hasFig('FigWav')
        hFigWav = obj.hFigs('FigWav');
        hFigWav.wait(0);
    end
end