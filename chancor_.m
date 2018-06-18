%--------------------------------------------------------------------------
% 9/22/17 JJJ: Created for SPARC
function mrCor = chancor_(mn, P)
    % mr = single(mn);
    mrCor = zeros(size(mn), 'single');
    nSites = numel(P.viSite2Chan);
    nSites_spk = P.maxSite*2+1-P.nSites_ref;
    miSites = P.miSites(2:nSites_spk,:);
    for iSite = 1:nSites
        viSites_ = miSites(:,iSite);
        vrDist_ = pdist2(P.mrSiteXY(iSite,:), P.mrSiteXY(viSites_,:));
        vrDist_ = 1 ./ (vrDist_ / min(vrDist_));
        vrCor_ = bsxfun(@times, single(mn(:,iSite)), single(mn(:,viSites_))) * vrDist_(:);
        mrCor(:,iSite) = conv(vrCor_, [1,2,3,2,1]/3, 'same');
        %     mrCor(:,iSrite) = vrCor_;
    end
    mrCor(mrCor<0)=0;
    mrCor = int16(-sqrt(mrCor));
    % figure;
    % subplot 121; imagesc(mrWav1(1:2000,:));
    % subplot 122; imagesc(((mrCor1(1:2000,:))));
end %func
