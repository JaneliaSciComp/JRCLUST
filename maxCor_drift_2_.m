%--------------------------------------------------------------------------
% 10/23/17 JJJ: find max correlation pair (combining drift and temporal shift)
function maxCor = maxCor_drift_2_(cmr1, cmr2, cviShift1, cviShift2)
    % assert(numel(cmr1) == numel(cmr2), 'maxCor_drift_: numel must be the same');
    nDrift = numel(cmr1);
    if nDrift == 1
        maxCor = max(xcorr2_mr_(cmr1{1}, cmr2{1}, cviShift1, cviShift2));
    else
        tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
        tr2 = cat(3, cmr2{:});
        nShift = numel(cviShift1);
        vrCor = zeros(1, nShift);
        for iShift = 1:nShift
            mr1_ = reshape(meanSubt_(tr1(cviShift1{iShift},:,:)), [], nDrift);
            mr2_ = reshape(meanSubt_(tr2(cviShift2{iShift},:,:)), [], nDrift);

            mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
            vrCor(iShift) = max(mr12_(:));
        end
        maxCor = max(vrCor);
    end
end
