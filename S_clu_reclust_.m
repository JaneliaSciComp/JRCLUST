%--------------------------------------------------------------------------
function S_clu = S_clu_reclust_(S_clu, S0, P);
    global trFet_spk
    vcMode_divide = 'amp'; % {'amp', 'density', 'fet'}

    trFet_spk0 = trFet_spk;
    nSites_fet = P.maxSite*2+1-P.nSites_ref;
    nFetPerSite = size(trFet_spk,1) / nSites_fet;
    switch vcMode_divide
        %     case 'nneigh'
        % %         % nearest neighbor averaging per same ecluster for feature enhancement
        %         trFet_spk = nneigh_ave_(S_clu, P, trFet_spk);
        %         P1 = setfield(P, 'nRepeat_fet', 1);
        %         S_clu = postCluster_(cluster_spacetime_(S0, P1), P);
        %         trFet_spk = trFet_spk0; % undo fet change (undo fet cleanup)

        %     case 'fet'
        % %         % recompute pca and
        %         vrSnr_clu = S_clu_snr_(S_clu);
        %         vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/nFetPerSite);
        %         vlRedo_spk = ismember(S_clu.viClu, find(vlRedo_clu));
        %         tnWav_spk = get_spkwav_(P, 0);
        %         trWav2_spk = single(permute(tnWav_spk(:,:,vlRedo_spk), [1,3,2]));
        %         trWav2_spk = spkwav_car_(trWav2_spk, P);
        %         [mrPv, vrD1] = tnWav2pv_(trWav2_spk, P);
        %         dimm1 = size(trWav2_spk);
        %         mrWav_spk1 = reshape(trWav2_spk, dimm1(1), []);
        %         trFet_spk_ = reshape(mrPv(:,1)' * mrWav_spk1, dimm1(2:3))';

        case 'amp'
        vrSnr_clu = S_clu_snr_(S_clu);
        try
            %         for iRepeat = (nFetPerSite-1):-1:1
            vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/2);
            vlRedo_spk = ismember(S_clu.viClu, find(vlRedo_clu));

            % reproject the feature
            %         nSpk_ = sum(vlRedo_spk);
            %         nFets_spk_ = ceil(size(trFet_spk,1)/2);
            %         trFet_spk_ = pca(reshape(trFet_spk(:,:,vlRedo_spk), size(trFet_spk,1), []), 'NumComponents', nFets_spk_);
            %         trFet_spk_ = permute(reshape(trFet_spk_, [size(trFet_spk,2), nSpk_, nFets_spk_]), [3,1,2]);
            %         trFet_spk = trFet_spk(1:nFets_spk_,:,:);
            %         trFet_spk(:,:,vlRedo_spk) = trFet_spk_;

            mlFet_ = false(nSites_fet, nFetPerSite);
            nSites_fet = ceil(nSites_fet*.5) %*.75
            mlFet_(1:nSites_fet, 1) = 1;
            %             mlFet_(1,:) = 1;
            %             mlFet_(:, 1) = 1;
            trFet_spk = trFet_spk0(find(mlFet_),:,:);

            %             trFet_spk = trFet_spk0(1:1*nSites_fet,:,:);
            S_clu_B = postCluster_(cluster_spacetime_(S0, P, vlRedo_spk), P);
            S_clu = S_clu_combine_(S_clu, S_clu_B, vlRedo_clu, vlRedo_spk);
        catch
            disperr_();
        end
        trFet_spk = trFet_spk0; %restore

        case 'density'
        vlRedo_clu = S_clu.vnSpk_clu > quantile(S_clu.vnSpk_clu, 1/2); %ilnear selection %2^(-iRepeat_clu+1)
        vlRedo_spk = ismember(S_clu.viClu, find(vlRedo_clu));
        S_clu_A = postCluster_(cluster_spacetime_(S0, P, ~vlRedo_spk), P);
        S_clu_B = postCluster_(cluster_spacetime_(S0, P, vlRedo_spk), P);
        S_clu.viClu(~vlRedo_spk) = S_clu_A.viClu;
        S_clu.viClu(vlRedo_spk) = S_clu_B.viClu + max(S_clu_A.viClu);

        otherwise
        return;
    end
end %func
