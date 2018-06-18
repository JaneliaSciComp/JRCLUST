%--------------------------------------------------------------------------
function trCov = tr2cov_(trWav)
    [nT, nSites, nClu] = size(trWav);
    trCov = zeros(nT, nT, nClu, 'like', trWav);
    trWav = meanSubt_(trWav);
    for iClu=1:nClu
        mrCov1 = trWav(:,:,iClu); %mean waveform covariance
        mrCov1 = mrCov1 * mrCov1';
        trCov(:,:,iClu) = mrCov1;
    end %for
end %func
