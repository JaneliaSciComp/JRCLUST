function resetFigTraces(hFigTraces, rawTraces, hCfg)
    %RESETFIGTRACES
    if hCfg.nSegmentsTraces > 1
        tlim1 = ([0, size(rawTraces, 2)] + hFigTraces.figData.windowBounds(1) - 1) / hCfg.sampleRate;
        tlim1 = round(tlim1*1000)/1000;

        hFigTraces.axApply('default', @axis, [tlim1, 0, hCfg.nSites + 1]);
    else
        hFigTraces.axApply('default', @axis, [hFigTraces.figData.windowBounds/hCfg.sampleRate, 0, hCfg.nSites + 1]);
    end
end
