%--------------------------------------------------------------------------
% 12/16/17 JJJ: Find overlapping spikes and set superthreshold sample points to zero in the overlapping region
function [spikeWaveforms_out, spikeWaveforms2_out] = cancel_overlap_spk_(spikeWaveforms, spikeWaveforms2, spikeTimes, spikeSites, spikeSecondarySites, siteThresholds, P)
    % Overlap detection. only return one stronger than other
    useGPU = isGpu_(spikeWaveforms);
    [spikeTimes, spikeWaveforms, spikeWaveforms2] = gather_(spikeTimes, spikeWaveforms, spikeWaveforms2);
    [viSpk_ol_spk, vnDelay_ol_spk, vnCount_ol_spk] = detect_overlap_spk_(spikeTimes, spikeSites, P);
    [spikeWaveforms_out, spikeWaveforms2_out] = deal(spikeWaveforms, spikeWaveforms2);
    % find spike index that are larger and fit and deploy
    viSpk_ol_a = find(viSpk_ol_spk>0); % later occuring
    [viSpk_ol_b, vnDelay_ol_b] = deal(viSpk_ol_spk(viSpk_ol_a), vnDelay_ol_spk(viSpk_ol_a)); % first occuring
    spikeTimes0 = int32(P.spkLim(1):P.spkLim(2));
    siteThresholds = gather_(-abs(siteThresholds(:))');
    % for each pair identify time range where threshold crossing occurs and set to zero
    % correct only first occuring (b)
    miSites = P.miSites;
    nSpk_ol = numel(viSpk_ol_a);
    nSpk = size(spikeWaveforms,2);
    for iSpk_ol = 1:nSpk_ol
        [iSpk_b, nDelay_b] = deal(viSpk_ol_b(iSpk_ol), vnDelay_ol_b(iSpk_ol));
        viSite_b = miSites(:,spikeSites(iSpk_b));
        mnWav_b = spikeWaveforms_out(nDelay_b+1:end,:,iSpk_b);
        mlWav_b = bsxfun(@le, mnWav_b, siteThresholds(viSite_b));
        mnWav_b(mlWav_b) = 0;
        spikeWaveforms_out(nDelay_b+1:end,:,iSpk_b) = mnWav_b;

        %     vnWav_b = median(spikeWaveforms(nDelay_b+1:end,:,iSpk_b), 2);
        %     spikeWaveforms_out(nDelay_b+1:end,:,iSpk_b) = repmat(vnWav_b, 1, nSpk);
        %     mnWav_b(nDelay_b+1:end,:) = repmat(vnWav_b, 1, nSpk);
        %     mlWav_b(1:nDelay_b,:) = 0; % safe time zone. remainders get cancelled
        %     mnWav_b(mlWav_b) = 0; % spike cancelled

        if ~isempty(spikeWaveforms2)
            viSite_b = miSites(:,spikeSecondarySites(iSpk_b));
            mnWav_b = spikeWaveforms2_out(nDelay_b+1:end,:,iSpk_b);
            mlWav_b = bsxfun(@le, mnWav_b, siteThresholds(viSite_b));
            mnWav_b(mlWav_b) = 0;
            spikeWaveforms2_out(nDelay_b+1:end,:,iSpk_b) = mnWav_b;
        end
    end %for
    % spikeWaveforms = gpuArray_(spikeWaveforms, useGPU);
    %     [iSite_a, iSite_b] = deal(spikeSites(iSpk_a), spikeSites(iSpk_b));
    %     [viSite_a, viSite_b] = deal(miSites(:,iSite_a), miSites(:,iSite_b));
    %     [viSite_ab, via_, vib_] = intersect(viSite_a, viSite_b);


    %     if isempty(viSite_ab), continue; end

    % find points under threshold
    %     [mnWav_a, mnWav_b] = deal(spikeWaveforms(:,via_,iSpk_a), spikeWaveforms(:,vib_,iSpk_b));
    %     vnThresh_ab = -siteThresholds(viSite_ab);
    %     [mlWav_a, mlWav_b] = deal(bsxfun(@lt, mnWav_a, vnThresh_ab), bsxfun(@lt, mnWav_b, vnThresh_ab));
    %

    %
    % correct both A and B by removing super-threshold points
    %     [mnWav_a, mnWav_b] = deal(spikeWaveforms(:,:,iSpk_a), spikeWaveforms(:,:,iSpk_b));
    %     [mlWav_a, mlWav_b] = deal(bsxfun(@lt, mnWav_a, -siteThresholds(viSite_a)), bsxfun(@lt, mnWav_b, -siteThresholds(viSite_b)));
    %     nDelay_b = vnDelay_ol_b(iSpk_ol);

    % set no overthreshold zone based on the delay, set it to half. only set superthreshold spikes to zero
end % function
