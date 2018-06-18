%--------------------------------------------------------------------------
% 10/23/17 JJJ: find max correlation pair (combining drift and temporal shift)
function maxCor = maxCor_drift_(cmr1, cmr2, cviShift1, cviShift2, fMode_cor)
    if nargin<5, fMode_cor=0; end %pearson corr
    assert_(numel(cmr1) == numel(cmr2), 'maxCor_drift_: numel must be the same');
    if numel(cmr1)==1
        maxCor = max(xcorr2_mr_(cmr1{1}, cmr2{1}, cviShift1, cviShift2));
    else
        tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
        tr2 = cat(3, cmr2{:});
        nDrift = numel(cmr1);
        nShift = numel(cviShift1);
        vrCor = zeros(1, nShift);
        for iShift = 1:nShift
            mr1_ = reshape(tr1(cviShift1{iShift},:,:), [], nDrift);
            mr2_ = reshape(tr2(cviShift2{iShift},:,:), [], nDrift);
            if fMode_cor == 0
                mr1_ = bsxfun(@minus, mr1_, sum(mr1_)/size(mr1_,1));
                mr2_ = bsxfun(@minus, mr2_, sum(mr2_)/size(mr2_,1));
                mr1_ = bsxfun(@rdivide, mr1_, sqrt(sum(mr1_.^2)));
                mr2_ = bsxfun(@rdivide, mr2_, sqrt(sum(mr2_.^2)));
                vrCor(iShift) = max(max(mr1_' * mr2_));
            else
                mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
                vrCor(iShift) = max(mr12_(:));
            end
        end
        maxCor = max(vrCor);

        %     tic
        %     mrCor = zeros(numel(cmr1), numel(cmr2));
        %     for i1=1:numel(cmr1)
        %         for i2=1:numel(cmr2)
        %             mrCor(i1,i2) = max(xcorr2_mr_(cmr1{i1}, cmr2{i2}, cviShift1, cviShift2));
        %         end
        %     end
        %     maxCor = max(mrCor(:));
        %     toc
    end
end
