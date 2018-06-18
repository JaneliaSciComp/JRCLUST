%--------------------------------------------------------------------------
function [mrWav_, mrWav_lo_, mrWav_hi_] = mean_wav_lo_hi_(tnWav_, viSpk_clu1, viSite_spk1, iSite_clu1, nSamples_max, P)

    vl_ = viSite_spk1 == iSite_clu1;
    viSpk_ = subsample_vr_(viSpk_clu1(vl_), nSamples_max);
    viSite_ = subsample_vr_(viSite_spk1(vl_), nSamples_max);
    mrWav_ = nanmean_int16_(tnWav_(:,:,viSpk_), 3, 1, iSite_clu1, viSite_, P);

    vl_lo = viSite_spk1 < iSite_clu1;
    if any(vl_lo)
        viSpk_lo_ = subsample_vr_(viSpk_clu1(vl_lo), nSamples_max);
        viSite_lo_ = subsample_vr_(viSite_spk1(vl_lo), nSamples_max);
        mrWav_lo_ = nanmean_int16_(tnWav_(:,:,viSpk_lo_), 3, 1, iSite_clu1, viSite_lo_, P);
    else
        mrWav_lo_ = mrWav_;
    end

    vl_hi = viSite_spk1 > iSite_clu1;
    if any(vl_hi)
        viSpk_hi_ = subsample_vr_(viSpk_clu1(vl_hi), nSamples_max);
        viSite_hi_ = subsample_vr_(viSite_spk1(vl_hi), nSamples_max);
        mrWav_hi_ = nanmean_int16_(tnWav_(:,:,viSpk_hi_), 3, 1, iSite_clu1, viSite_hi_, P);
    else
        mrWav_hi_ = mrWav_;
    end
end %func
