%--------------------------------------------------------------------------
function [mrPv, vrD1] = tnWav2pv_(tr, P)
    %tr: nSamples x nSpikes x nChans
    MAX_SAMPLE = 10000;

    if nargin<2, P = get0_('P'); end
    % if nargin<3, viSites_ref = []; end
    % if isempty(tr)
    %     nSpk = size(tnWav_spk,3);
    %     viSpk_sub = subsample_vr_(1:nSpk, MAX_SAMPLE);
    %     tr = permute(tnWav_spk(:,:,viSpk_sub), [1 3 2]);
    %     tr = single(tr);
    %     tr = spkwav_car_(tr, P);
    %     mrSpkWav1 = tr(:,:,1);
    % else
    viSpk_sub = subsample_vr_(1:size(tr,2), MAX_SAMPLE);
    mrSpkWav1 = tr(:,viSpk_sub, 1);
    % end
    % tr = single(tr);
    % if ~isempty(viSites_ref), tr = spkwav_car_(tr, viSites_ref); end

    % mrCov = meanSubt_(mrSpkWav1);
    mrCov = mrSpkWav1 * mrSpkWav1';
    [mrPv1, vrD1] = eig(mrCov);
    mrPv1 = jrclust.utils.zscore(fliplr(mrPv1)); % sort largest first
    vrD1 = flipud(diag(vrD1));

    % spike center should be negative
    iMid = 1-P.spkLim(1);
    vrSign = (mrPv1(iMid,:) < 0) * 2 - 1; %1 or -1 depending on the sign
    mrPv = bsxfun(@times, mrPv1, vrSign);
end
