%--------------------------------------------------------------------------
% 10/23/17 JJJ: find 1 - max(normalized distance), combining drift and temporal shift
function minDist = minDist_drift_(cmr1, cmr2, cviShift1, cviShift2)
    assert_(numel(cmr1) == numel(cmr2), 'minDist_drift_: numel must be the same');
    if numel(cmr1)==1
        a_ = reshape(cmr1{1},[],1); % TW?
        b_ = reshape(cmr2{1},[],1); % TW?

        minDist = 1 - sum(abs(a_-b_))/max([sum(abs(a_)); sum(abs(b_))]); % TW?
    else
        tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
        tr2 = cat(3, cmr2{:});

        % begin TW block
        c_ = zeros(size(tr1, 3), 1);
        for iDrift=1:size(tr1,3)
            a_ = reshape(tr1(:,:,iDrift),[],1);
            b_ = reshape(tr2(:,:,iDrift),[],1);
            c_(iDrift) = 1 - sum(abs(a_-b_))/max([sum(abs(a_)); sum(abs(b_))]);
        end
            
        minDist = max(c_);
        % end TW block
    end
end
