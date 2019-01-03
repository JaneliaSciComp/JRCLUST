function resetFigTraces(hFigTraces, rawTraces, hCfg)
    %RESETFIGTRACES
    nSites = numel(hCfg.siteMap);

    if hCfg.nTime_traces > 1
        tlim1 = ([0, size(rawTraces, 2)] + hFigTraces.figData.windowBounds(1) - 1) / hCfg.sampleRate;
        tlim1 = round(tlim1*1000)/1000;

        hFigTraces.axis([tlim1, 0, nSites + 1]);
    else
        hFigTraces.axis([hFigTraces.figData.windowBounds/hCfg.sampleRate, 0, nSites + 1]);
    end
end
