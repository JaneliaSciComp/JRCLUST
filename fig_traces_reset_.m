%--------------------------------------------------------------------------
function fig_traces_reset_(S_fig)
    global mnWav1

    if nargin<1, [hFig, S_fig] = getCachedFig('Fig_traces');  end
    % axis(S_fig.hAx, [S_fig.nlim_bin / P.sampleRateHz, 0, nSites+1]);
    P = get0_('P');
    nTime_traces = get_(P, 'nTime_traces');
    if nTime_traces > 1
        tlim1 = ([0, size(mnWav1,1)] + S_fig.nlim_bin(1) - 1) / P.sampleRateHz;
        tlim1 = round(tlim1*1000)/1000;
        axis_(S_fig.hAx, [tlim1, 0, numel(P.chanMap)+1]);
    else
        axis_(S_fig.hAx, [S_fig.nlim_bin / P.sampleRateHz, 0, numel(P.chanMap)+1]);
    end
end % function
