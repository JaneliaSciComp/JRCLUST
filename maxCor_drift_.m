%--------------------------------------------------------------------------
% 10/23/17 JJJ: find max correlation pair (combining drift and temporal shift)
function maxCor = maxCor_drift_(cmr1, cmr2, cviShift1, cviShift2, fMode_cor)
    if nargin<5, fMode_cor=0; end %pearson corr
    dialogAssert(numel(cmr1) == numel(cmr2), 'maxCor_drift_: numel must be the same');
    if numel(cmr1)==1
        maxCor = max(xcorr2_mr_(cmr1{1}, cmr2{1}, cviShift1, cviShift2));
        
        a_=reshape(cmr1{1},[],1); % TW?
        b_=reshape(cmr2{1},[],1); % TW?

        % m_=max([var(a_) var(b_)]); % TW?
        % tmp_=cov(a_,b_)/m_; % TW?
        % c_(counter) = tmp_(1,2); % TW?
        maxCor= 1-sum(abs(a_-b_))/max([sum(abs(a_));sum(abs(b_))]); % TW?
    else
        tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
        tr2 = cat(3, cmr2{:});
        % nDrift = numel(cmr1);
        % nShift = numel(cviShift1);
        % vrCor = zeros(1, nShift);
        % for iShift = 1:nShift
        %     mr1_ = reshape(tr1(cviShift1{iShift},:,:), [], nDrift);
        %     mr2_ = reshape(tr2(cviShift2{iShift},:,:), [], nDrift);
        %     if fMode_cor == 0
        %         mr1_ = bsxfun(@minus, mr1_, sum(mr1_)/size(mr1_,1));
        %         mr2_ = bsxfun(@minus, mr2_, sum(mr2_)/size(mr2_,1));
        %         mr1_ = bsxfun(@rdivide, mr1_, sqrt(sum(mr1_.^2)));
        %         mr2_ = bsxfun(@rdivide, mr2_, sqrt(sum(mr2_.^2)));
        %         vrCor(iShift) = max(max(mr1_' * mr2_));
        %     else
        %         mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
        %         vrCor(iShift) = max(mr12_(:));
        %     end
        % end
        % maxCor = max(vrCor);

        % begin TW block
        c_ = zeros(size(tr1, 3), 1);
        for counter=1:size(tr1,3)
            a_=reshape(tr1(:,:,counter),[],1);
            b_=reshape(tr2(:,:,counter),[],1);
        
        % m_=max([var(a_) var(b_)]);
        % tmp_=cov(a_,b_)/m_;
        % c_(counter) = tmp_(1,2);  
        c_(counter) = -sum(abs(a-b))/max(abs([a_ b_]));
        end
            
        maxCor=max(c_);
        % end TW block

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
