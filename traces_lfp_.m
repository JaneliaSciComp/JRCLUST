%--------------------------------------------------------------------------
% Show LFP traces
function traces_lfp_(P)
    P.sampleRateHz = P.sampleRateHz_lfp;
    P.vcFile = P.lfpFile;
    P.multiFilenames = {};
    P.tlim(2) = P.tlim(1) + diff(P.tlim) * P.nSkip_lfp;
    P.headerOffset = 0;
    traces_(P, 0, [], 1); % don't show spikes, don't set P to S0
end %func
