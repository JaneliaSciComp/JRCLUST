%--------------------------------------------------------------------------
function S_clu = S_clu_reclust_(S_clu, S0, P);
    global spikeFeatures

    if ~all(size(spikeFeatures) ~= S0.featureDims)
        spikeFeatures = getSpikeFeatures(P);
    end

    vcMode_divide = 'amp'; % {'amp', 'density', 'fet'}

    spikeFeatures0 = spikeFeatures;
    nSites_fet = P.maxSite*2+1-P.nSites_ref;
    nFetPerSite = size(spikeFeatures,1) / nSites_fet;
    switch vcMode_divide
        %     case 'nneigh'
        % %         % nearest neighbor averaging per same ecluster for feature enhancement
        %         spikeFeatures = nneigh_ave_(S_clu, P, spikeFeatures);
        %         P1 = setfield(P, 'nRepeat_fet', 1);
        %         S_clu = postCluster_(cluster_spacetime_(S0, P1), P);
        %         spikeFeatures = spikeFeatures0; % undo fet change (undo fet cleanup)

        %     case 'fet'
        % %         % recompute pca and
        %         vrSnr_clu = S_clu_snr_(S_clu);
        %         vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/nFetPerSite);
        %         vlRedo_spk = ismember(S_clu.spikeClusters, find(vlRedo_clu));
        %         spikeWaveforms = getSpikeWaveforms(P, 0);
        %         trWav2_spk = single(permute(spikeWaveforms(:,:,vlRedo_spk), [1,3,2]));
        %         trWav2_spk = spkwav_car_(trWav2_spk, P);
        %         [mrPv, vrD1] = tnWav2pv_(trWav2_spk, P);
        %         dimm1 = size(trWav2_spk);
        %         mrWav_spk1 = reshape(trWav2_spk, dimm1(1), []);
        %         spikeFeatures_ = reshape(mrPv(:,1)' * mrWav_spk1, dimm1(2:3))';

        case 'amp'
        vrSnr_clu = S_clu_snr_(S_clu);
        try
            %         for iRepeat = (nFetPerSite-1):-1:1
            vlRedo_clu = vrSnr_clu < quantile(vrSnr_clu, 1/2);
            vlRedo_spk = ismember(S_clu.spikeClusters, find(vlRedo_clu));

            % reproject the feature
            %         nSpk_ = sum(vlRedo_spk);
            %         nFets_spk_ = ceil(size(spikeFeatures,1)/2);
            %         spikeFeatures_ = pca(reshape(spikeFeatures(:,:,vlRedo_spk), size(spikeFeatures,1), []), 'NumComponents', nFets_spk_);
            %         spikeFeatures_ = permute(reshape(spikeFeatures_, [size(spikeFeatures,2), nSpk_, nFets_spk_]), [3,1,2]);
            %         spikeFeatures = spikeFeatures(1:nFets_spk_,:,:);
            %         spikeFeatures(:,:,vlRedo_spk) = spikeFeatures_;

            mlFet_ = false(nSites_fet, nFetPerSite);
            nSites_fet = ceil(nSites_fet*.5) %*.75
            mlFet_(1:nSites_fet, 1) = 1;
            %             mlFet_(1,:) = 1;
            %             mlFet_(:, 1) = 1;
            spikeFeatures = spikeFeatures0(find(mlFet_),:,:);

            %             spikeFeatures = spikeFeatures0(1:1*nSites_fet,:,:);
            S_clu_B = postCluster_(cluster_spacetime_(S0, P, vlRedo_spk), P);
            S_clu = S_clu_combine_(S_clu, S_clu_B, vlRedo_clu, vlRedo_spk);
        catch
            disperr_();
        end
        spikeFeatures = spikeFeatures0; %restore

        case 'density'
        vlRedo_clu = S_clu.nSpikesPerCluster > quantile(S_clu.nSpikesPerCluster, 1/2); %ilnear selection %2^(-iRepeat_clu+1)
        vlRedo_spk = ismember(S_clu.spikeClusters, find(vlRedo_clu));
        S_clu_A = postCluster_(cluster_spacetime_(S0, P, ~vlRedo_spk), P);
        S_clu_B = postCluster_(cluster_spacetime_(S0, P, vlRedo_spk), P);
        S_clu.spikeClusters(~vlRedo_spk) = S_clu_A.viClu;
        S_clu.spikeClusters(vlRedo_spk) = S_clu_B.viClu + max(S_clu_A.viClu);

        otherwise
        return;
    end
end % function
