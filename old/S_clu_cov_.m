%--------------------------------------------------------------------------
function mrDist_clu = S_clu_cov_(S_clu, P)
    % merge clusters based on the spike waveforms
    % nPc = 3;
    % trCov = tr2cov_(S_clu.trWav_spk_clu);
    % trWav_clu = meanSubt_(S_clu.trWav_spk_clu);
    trWav_clu = meanSubt_(S_clu.trWav_raw_clu); %uses unfiltered waveform
    % pairwise similarity computation
    % maxSite = P.maxSite_merge;
    maxSite = ceil(P.maxSite/2);
    % maxSite = P.maxSite;
    mrDist_clu = nan(S_clu.nClu);
    vrVar_clu = zeros(S_clu.nClu, 1, 'single');
    % miSites = P.miSites(1:P.maxSite_merge*2+1, :);
    for iClu2 = 1:S_clu.nClu
        iSite2 = S_clu.viSite_clu(iClu2);
        %     viSite2 = miSites(:,iSite2);
        viClu1 = find(abs(S_clu.viSite_clu - iSite2) <= maxSite);
        viClu1(viClu1 == iClu2) = [];
        viClu1 = viClu1(:)';
        %     mrCov2 = eigvec_(trCov(:,:,iClu2), nPc);
        mrWav_clu2 = trWav_clu(:,:,iClu2);
        [~,~,var2] = pca(mrWav_clu2); var2 = var2(1)/sum(var2);
        vrVar_clu(iClu2) = var2;
        for iClu1 = viClu1
            %         mrCov1 = eigvec_(trCov(:,:,iClu1), nPc);
            mrWav_clu1 = trWav_clu(:,:,iClu1);
            %         mrDist_clu(iClu1, iClu2) = mean(abs(mean(mrCov1 .* mrCov2))); %cov_eig_(mrCov1, mrCov2, nPc);
            mrDist_clu(iClu1, iClu2) = cluWav_dist_(mrWav_clu1, mrWav_clu2);
        end
    end
end %func
