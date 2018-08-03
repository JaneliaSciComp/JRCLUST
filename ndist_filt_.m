%--------------------------------------------------------------------------
function mnWav2 = ndist_filt_(mnWav2, ndist_filt)

    vnFilt_ = gpuArray_(ones(ndist_filt,1,'single'), isGpu_(mnWav2));
    mnWav_ = mnWav2(1+ndist_filt:end,:) - mnWav2(1:end-ndist_filt,:);
    [n1, nChans] = deal(round((ndist_filt-1)/2) , size(mnWav_,2));
    n2 = ndist_filt - n1;
    mnWav_ = [zeros([n1, nChans], 'like', mnWav_); mnWav_; zeros([n2, nChans], 'like', mnWav_)];
    dataType_ = class_(mnWav2);
    for iChan = 1:nChans
        vn_ = cast(sqrt(conv(single(mnWav_(:,iChan)).^2, vnFilt_, 'same')), dataType_);
        vn_ = median(vn_(1:10:end)) - vn_;
        mnWav2(:,iChan) = vn_;
    end
end %func
