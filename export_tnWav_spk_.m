%--------------------------------------------------------------------------
function export_spikeWaveforms_(h,e)
    % Export all spike waveforms from selected cluster

    S0 = get(0, 'UserData');
    P = S0.P;
    spikeWaveforms = get_spkwav_(P, 0);
    spikeTraces = get_spkwav_(P, 1);
    S_clu = S0.S_clu;
    iClu1 = S0.iCluCopy;
    viSpk1 = S_clu.cviSpk_clu{iClu1};
    nSpk1 = numel(viSpk1);
    nSites = numel(P.chanMap);
    waveformDims1 = size(spikeWaveforms);    waveformDims1(2) = nSites;
    traceDims1 = size(spikeTraces);    traceDims1(2) = nSites;
    spikeWaveforms1 = zeros(waveformDims1, 'like', spikeWaveforms);
    spikeTraces1 = zeros(traceDims1, 'like', spikeTraces);
    miSites_spk1 = P.miSites(:, S0.spikeSites);
    for iSpk = 1:nSpk1
        spikeWaveforms1(:, miSites_spk1(:,iSpk), iSpk) = spikeWaveforms(:,:,iSpk);
        spikeTraces1(:, miSites_spk1(:,iSpk), iSpk) = spikeTraces(:,:,iSpk);
    end
    eval(sprintf('spikeWaveforms_clu%d = spikeWaveforms1;', iClu1));
    eval(sprintf('assignWorkspace_(spikeWaveforms_clu%d);', iClu1));
    eval(sprintf('spikeTraces_clu%d = spikeTraces1;', iClu1));
    eval(sprintf('assignWorkspace_(spikeTraces_clu%d);', iClu1));
end %func
