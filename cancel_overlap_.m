%--------------------------------------------------------------------------
% 2017/12/15 JJJ: correct overlapping spikes
% Question: how to best fit? using global spike pca (pick) or cluster-private pca
function trFet_spk = cancel_overlap_(cviSpk_o_1, cviSpk_o_12, cviDelay1, S0, S_clu, P)
    global tnWav_spk trFet_spk tnWav_raw

    viSites_ref = ceil(size(tnWav_spk,2)/2):size(tnWav_spk,2);
    nPc_fit = min(get_set_(P, 'nPc_fit_overlap', 5), size(S0.mrPv_global,2));
    mrPv = S0.mrPv_global(:,1:nPc_fit) / sqrt(size(S0.mrPv_global,1));
    % fit and cancel overlap using delay and mean clu waveform
    for iClu = 1:S_clu.nClusters
        [viSpk_o_1, viSpk_o_12, viDelay1] = deal(cviSpk_o_1{iClu}, cviSpk_o_12{iClu}, cviDelay1{iClu});
        if isempty(viSpk_o_1), continue; end
        [mrWav_clu1, iSite1] = deal(S_clu.trWav_spk_clu(:,:,iClu), S_clu.clusterSites(iClu));
        [clusterSites1, viSite1, viSite12] = deal(P.miSites(:, iSite1), S0.spikeSites(viSpk_o_1), S0.spikeSites(viSpk_o_12));
        trWav1 = trWav_car_(tnWav_spk(:,:,viSpk_o_1), P);
        trWav12 = trWav_car_(tnWav_spk(:,:,viSpk_o_12), P);
        %     [tnWav1, tnWav12] = deal(tnWav_spk(:,:,viSpk_o_1), tnWav_spk(:,:,viSpk_o_12));
        %     mrWav_clu1_car = bsxfun(@minus, mrWav_clu1, mean(mrWav_clu1(:,viSites_ref),2));
        %     mrWav_clu1_shift = vr_shift_(mrWav_clu1_car(:,1), -1:.5:1);
        % fft and use power to fit waveforms using the primary sites only. preserve phase and change power

        % calculate private pc. using non-colliding spikes (?). generate basis for the cluster
        %     viSpk_clu1 = S_clu.cviSpk_clu{iClu};
        %     viSpk_clu1 = subsample_vr_(viSpk_clu1(S0.spikeSites(viSpk_clu1) == iSite1), 5e3);
        %     mr_ = reshape(meanSubt_(tnWav_spk(:,:,viSpk_clu1)), [], numel(viSpk_clu1));
        %     [~,mrPv_clu1] = pca(mr_, 'NumComponents', 3); %subsample if must

        % fit
        mrPv = mrPv_clu1;

        %     trWav1 = single(tnWav1);
        iSite_fit = 1:14;
        vrData = toVec(meanSubt_(trWav1(:,iSite_fit,1)));
        a = mrPv' * vrData;  % coefficient
        mrWav1_fit = mrPv * a;
        mrWav_err1 = trWav1(:,iSite_fit,1) - mrWav1_fit;
        myfig; plot(meanSubt_(trWav1(:,iSite_fit,1)), 'k'); plot(mrWav_err1, 'r');
        myfig; plot(mrWav1_fit, 'g'); plot(mrWav_err1, 'r');
        %     tnWav_raw1 = tnWav_raw(:,:,viSpk_o_1);
        %     tnWav_raw12 = tnWav_raw(:,:,viSpk_o_12);
        % fit using center channel only
        %     trWav_car1 = trWav_car_(trWav1, P);
        %     trWav_car12 = trWav_car_(trWav12, P);

    end
end %func
