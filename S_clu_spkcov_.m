%--------------------------------------------------------------------------
function mrDist_clu = S_clu_spkcov_(S_clu, P)
    % compute covariance from each spikes belonging to clusters
    MAX_SAMPLE = 1000;
    nPc = 3;
    viSite_spk = get0_('viSite_spk');
    maxSite = ceil(P.maxSite/2);
    mrDist_clu = nan(S_clu.nClu);
    tnWav_raw = get_spkwav_(P, 1);

    % compute a cell of mrPv
    [cmrPv_spk_clu, cmrPv_raw_clu] = deal(cell(S_clu.nClu, 1));
    for iClu = 1:S_clu.nClu
        viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        viSpk_clu1 = viSpk_clu1(viSite_spk(viSpk_clu1) == S_clu.viSite_clu(iClu));
        viSpk_clu1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
        %     cmrPv_spk_clu{iClu} = tn2pca_spk_(tnWav_spk(:,:,viSpk_clu1), nPc);
        cmrPv_raw_clu{iClu} = tn2pca_spk_(tnWav_raw(:,:,viSpk_clu1), nPc);
    end %for

    for iClu2 = 1:S_clu.nClu
        iSite2 = S_clu.viSite_clu(iClu2);
        viClu1 = find(abs(S_clu.viSite_clu - iSite2) <= maxSite);
        viClu1(viClu1 == iClu2) = [];
        viClu1 = viClu1(:)';
        %     mrPv_clu2 = cmrPv_clu_raw{iClu2};
        for iClu1 = viClu1
            vrCorr_raw = diag(corr_(cmrPv_raw_clu{iClu2}, cmrPv_raw_clu{iClu1}));
            %         vrCorr_spk = diag(corr_(cmrPv_spk_clu{iClu2}, cmrPv_spk_clu{iClu1}));
            %         mrDist_clu(iClu1, iClu2) = mean(abs([vrCorr_raw; vrCorr_spk]));
            %         mrDist_clu(iClu1, iClu2) = mean(abs([vrCorr_raw; vrCorr_spk]));
            %         mrPv_clu1 = cmrPv_clu_raw{iClu1};
            mrDist_clu(iClu1, iClu2) = mean(abs(vrCorr_raw));
        end
    end %for
end %func
