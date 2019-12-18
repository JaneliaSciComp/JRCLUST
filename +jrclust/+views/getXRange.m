function xRange = getXRange(iCluster, nSamples, hCfg)
    %GETXRANGE Get the x-values of iCluster waveforms for main view
    if hCfg.showRaw
        evtWindowSamp = hCfg.evtWindowRawSamp;
    else
        evtWindowSamp = hCfg.evtWindowSamp;
    end

    xOffset = iCluster - 1 + evtWindowSamp(2)/nSamples;

    xRange = xOffset + (1:nSamples)/nSamples;
    xRange([1, end]) = nan;
    xRange = single(xRange(:));
end
