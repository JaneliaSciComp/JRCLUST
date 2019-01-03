%--------------------------------------------------------------------------
% Show LFP traces
function traces_lfp_(P)
    P.sRateHz = P.sRateHz_lfp;
    P.vcFile = P.vcFile_lfp;
    P.csFile_merge = {};
    P.tlim(2) = P.tlim(1) + diff(P.tlim) * P.nSkip_lfp;
    P.header_offset = 0;
    traces_(P, [], 1); % don't show spikes, don't set P to S0
end %func
