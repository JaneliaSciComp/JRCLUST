%--------------------------------------------------------------------------
function selfcorr = S_clu_self_corr_(S_clu, iClu1, S0)
    % plot top half vs bottom half correlation. sum of vpp
    if nargin<2, iClu1 = []; end
    if nargin<3, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    [spikeSites, P] = deal(S0.spikeSites, S0.P);
    spikeTraces = get_spkwav_(P, getOr(P, 'fWavRaw_merge', 1));

    if isempty(iClu1)
        fprintf('Computing self correlation\n\t'); t1=tic;
        selfcorr = zeros(1, S_clu.nClusters);
        for iClu=1:S_clu.nClusters
            selfcorr(iClu) = S_clu_self_corr__(S_clu, spikeTraces, iClu, spikeSites);
            fprintf('.');
        end
        fprintf('\n\ttook %0.1fs\n', toc(t1));
    else
        selfcorr = S_clu_self_corr__(S_clu, spikeTraces, iClu1, spikeSites);
    end
end %func