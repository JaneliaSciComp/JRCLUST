%--------------------------------------------------------------------------
function [mrPc1, mrPc2] = pca_tr_(tn)
    % returns first principal component across sites
    % persistent tn_
    % global tnWav_spk
    %
    % if nargin<1, tn = tn_; end
    % if isempty(tn)
    %     tn_ = tnWav_spk;
    %     tn = tn_;
    % end

    % mrPv = zeros(size(tr,1), size(tr,3), 'single');
    mrPc1 = zeros(size(tn,2), size(tn,3), 'single');
    mrPc2 = ifeq_(nargout>1, mrPc1, []);

    % tr = single(tr);
    % n = size(tn,2);
    tr = meanSubt_tr_(single(tn));
    nSpk = size(tn,3);
    % tic
    if isempty(mrPc2)
        parfor iSpk = 1:nSpk
            mr1 = tr(:,:,iSpk);
            [V,D] = eig(mr1*mr1');
            D = sqrt(diag(D));
            mrPc1(:,iSpk) = V(:,end)' * mr1 / D(end); %V(:,end)' * mr1;
        end
    else
        parfor iSpk = 1:nSpk
            mr1 = tr(:,:,iSpk);
            [V,D] = eig(mr1*mr1');
            D = sqrt(diag(D));
            mrPc1(:,iSpk) = V(:,end)' * mr1 / D(end); %V(:,end)' * mr1;
            mrPc2(:,iSpk) = V(:,end-1)' * mr1 / D(end-1); %V(:,end)' * mr1;
        end
    end
    % toc
end %func
