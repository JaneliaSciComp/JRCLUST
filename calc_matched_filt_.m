%--------------------------------------------------------------------------
function [vrFilt_spk, vrVaf, nShift_post] = calc_matched_filt_(mnWav1, P) %detect primary
    % generate a matched filter Kernel
    % determine the treshold
    % 6/29/17 JJJ: Spike waveform matched fitler determination
    % 6/30/17 JJJ: Range optimization

    vnThresh_site = gather_(int16(mr2rms_(mnWav1, 1e5) * P.qqFactor));
    [viTime_spk, vnAmp_spk, viSite_spk] = detect_spikes_(mnWav1, vnThresh_site, [], P);

    % extract wave forms
    nSpks = numel(viSite_spk);
    nSites = numel(P.chanMap);
    %spkLim = [-1, 1] * round(mean(abs(P.spkLim))); %balanced
    % spkLim = [-1, 1] * round(max(abs(P.spkLim))); %balanced
    % spkLim = [-1, 1] * round(min(abs(P.spkLim))); %balanced
    % spkLim = spkLim + [1,0];
    spkLim = [-1,1] * round(mean(abs(P.spkLim)));

    mnWav_spk = zeros(diff(spkLim) + 1, nSpks, 'int16');
    for iSite = 1:nSites
        viiSpk11 = find(viSite_spk == iSite);
        if isempty(viiSpk11), continue; end
        viTime_spk11 = viTime_spk(viiSpk11); %already sorted by time
        mnWav_spk(:,viiSpk11) = gather_(vr2mr3_(mnWav1(:,iSite), viTime_spk11, spkLim));
    end

    [vrFilt_spk, ~, vrVaf] = pca(single(mnWav_spk'), 'NumComponents', 1, 'Centered', 0);
    vrFilt_spk = flipud(vrFilt_spk(:));
    if abs(min(vrFilt_spk)) > abs(max(vrFilt_spk)), vrFilt_spk = -vrFilt_spk; end
    % vrFilt_spk(1) = []; % start from 1 correction
    % vrFilt_spk = vrFilt_spk - mean(vrFilt_spk);

    % vrFilt_spk = vrFilt_spk / (vrFilt_spk.'*vrFilt_spk);
    vrVaf = cumsum(vrVaf);
    vrVaf = vrVaf / vrVaf(end);
    [~,nShift_post] = max(vrFilt_spk);
    nShift_post = round(numel(vrFilt_spk)/2 - nShift_post);
    % if nShift_post < 0
    %     vrFilt_spk = vrFilt_spk(1-nShift_post:end);
    % else
    %     vrFilt_spk = vrFilt_spk(1:end-nShift_post);
    % end
    % nShift_post = 0;
end %func
