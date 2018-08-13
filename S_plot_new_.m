%--------------------------------------------------------------------------
function S_plot = S_plot_new_(S0)
    % S_plot contains quantities to be plotted
    % Copied from jrclust.m quality_metric_

    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P;
    spikeWaveforms = get_spkwav_(P, 0);

    vrVrms_site = single(S0.vrThresh_site(:)) / P.qqFactor;
    vrSnr_evt = single(abs(S0.vrAmp_spk(:))) ./ vrVrms_site(S0.spikeSites(:));
    t_dur = double(max(S0.spikeTimes) - min(S0.spikeTimes)) / P.sRateHz;
    vrRate_site = cellfun(@numel, S0.cviSpk_site)' / t_dur;
    % nSites = numel(S0.cviSpk_site);

    % calc # spikes exceeding detection threshold
    vnSite_evt = zeros(size(S0.spikeTimes), 'int16');
    for iSite = 1:numel(S0.cviSpk_site)
        viSpk_site1 = S0.cviSpk_site{iSite};
        mrMin_site1 = squeeze_(min(spikeWaveforms(:,:,viSpk_site1)));
        vrThresh_site1 = -abs(S0.vrThresh_site(P.miSites(:, iSite)));
        vnSite_evt(viSpk_site1) = sum(bsxfun(@lt, mrMin_site1, vrThresh_site1(:)));
    end

    % cluster meta analysis (cluster of clusters)


    % Compute cluster stats
    mrMin_clu = uV2bit_(squeeze_(min(S0.S_clu.trWav_spk_clu)));
    vrSnr_clu = abs(mrMin_clu(1,:))' ./ vrVrms_site(S0.S_clu.clusterSites);
    vrRate_clu = cellfun(@numel, S0.S_clu.spikesByCluster)' / t_dur;
    mrThresh_clu = -abs(S0.vrThresh_site(P.miSites(:,S0.S_clu.clusterSites)));
    vnSite_clu = sum(mrMin_clu < mrThresh_clu)';

    S_plot = makeStruct_(vrVrms_site, vrRate_site, t_dur, P, ...
    vrSnr_evt, vnSite_evt, vrSnr_clu, vrRate_clu, vnSite_clu);
end %func
