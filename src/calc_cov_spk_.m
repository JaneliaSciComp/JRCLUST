%--------------------------------------------------------------------------
function [mrVpp1, mrVpp2] = calc_cov_spk_(viSpk1, viSites1)

    [viSite_spk, P] = get0_('viSite_spk', 'P');
    tnWav_spk = get_spkwav_(P, 0); % get filtered waveform

    nSpk1 = numel(viSpk1);
    viSites_spk1 = viSite_spk(viSpk1);
    tnWav_spk1 = gpuArray_(tnWav_spk(:,:,viSpk1), P.fGpu);
    nSites_spk = 1 + P.maxSite * 2;
    [mrVpp1_, mrVpp2_] = trWav2fet_(tnWav_spk1, P, nSites_spk);
    [mrVpp1_, mrVpp2_] = multifun_(@(x)gather_(abs(x)), mrVpp1_, mrVpp2_);

    % re-project to common basis
    viSites_spk_unique = unique(viSites_spk1);
    [mrVpp1, mrVpp2] = deal(zeros([numel(viSites1), nSpk1], 'like', mrVpp1_));
    for iSite1 = 1:numel(viSites_spk_unique) %only care about the first site
        iSite11 = viSites_spk_unique(iSite1); %center sites group
        viSpk11 = find(viSites_spk1 == iSite11); %dangerous error
        viSites11 = P.miSites(:, iSite11);
        [vlA11, viiB11] = ismember(viSites11, viSites1);
        mrVpp1(viiB11(vlA11),viSpk11) = mrVpp1_(vlA11,viSpk11);
        mrVpp2(viiB11(vlA11),viSpk11) = mrVpp2_(vlA11,viSpk11);
    end
end %func
