%--------------------------------------------------------------------------
% function [S_clu, nClu_merged] = S_clu_pv_merge_(S_clu, P) %update mrWavCor when you merge
%
% MAD_THRESH = -4;
% nClu_merged = 0;
% % [viSite_spk] = get0_('viSite_spk');
% % mrWavCor = S_clu.mrWavCor;
% nClu = S_clu.nClu;
% fprintf('S_clu_pv_merge_\n');
%
% % Identify clusters to remove, update and same (no change), disjoint sets
% [vrMinDist_logz_clu, viMinDist_clu] = S_clu_pca_dist_(S_clu);
% %vi_clu1 = find(vrMinDist_logz_clu < MAD_THRESH);
% vi_clu1 = find(vrMinDist_logz_clu < -max(vrMinDist_logz_clu));
% if isempty(vi_clu1), return; end
%
% vi_clu2 = viMinDist_clu(vi_clu1);
% viMap_clu = 1:nClu;
% viMap_clu(vi_clu1) = vi_clu2;
% viClu_same = setdiff(1:nClu, union(vi_clu1, vi_clu2));
% viClu_remove = setdiff(1:nClu, viMap_clu);
% viClu_update = setdiff(setdiff(1:nClu, viClu_same), viClu_remove);
% % viClu_update = setdiff(1:nClu, viClu_same);
%
% % update cluster number
% try S_clu.icl(viClu_remove) = []; catch, end
% S_clu = S_clu_map_index_(S_clu, viMap_clu); %index mapped
% P.fVerbose = 0;
% S_clu = S_clu_refrac_(S_clu, P); % remove refrac spikes
%
% % update cluster waveforms and distance
% S_clu = S_clu_wav_(S_clu, viClu_update); %update cluster waveforms
% S_clu.mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update);
% % S_clu = S_clu_refresh_(S_clu); % remove empty and remap
% S_clu = S_clu_remove_empty_(S_clu);
%
% nClu_merged = nClu - S_clu.nClu;
% fprintf('\n\tnClu: %d->%d (%d merged)\n', nClu, S_clu.nClu, nClu_merged);
% end %func


%--------------------------------------------------------------------------
% function [vrMinDist_logz_clu, viMinDist_clu] = S_clu_pca_dist_(S_clu)
% global tnWav_raw tnWav_spk
%
% MAX_REAL_DIST = 50;
% MAX_SAMPLE = 2000;
% fUseMean = 1; %use median instead
% nPc = 2;
% fUseRaw = 1;
% fUsePvCorr = 1;
% nShift = 6;
% fUseSd = 1;
%
% P = S_clu.P;
%
% trWav_clu = ifeq_(fUseRaw, S_clu.trWav_raw_clu, S_clu.trWav_spk_clu);
% if ~fUseMean
%     viSite_spk = get0_('viSite_spk');
% end
% nClu = S_clu.nClu;
% nSamples = size(trWav_clu,1);
% % mrPv1_clu = zeros(nSamples, nClu);
% % mrPv1_clu = zeros(nSamples, nClu);
% nDelay = 3;
% [mrPv1_clu, mrPv2_clu, mrPv3_clu] = deal(zeros(size(trWav_clu,1), nClu));
% for iClu=1:nClu
% %     [~, mrPv1_clu(:,iClu)] = pca(trWav_clu(:,:,iClu), 'NumComponents', 1);
%     if fUseMean
%         mrWav_clu1 = trWav_clu(:,:,iClu);
%     else
%         viSpk_clu1 = S_clu.cviSpk_clu{iClu};
%         viSpk_clu1 = viSpk_clu1(viSite_spk(viSpk_clu1) == S_clu.viSite_clu(iClu));
%         viSpk_clu1 = subsample_vr_(viSpk_clu1, MAX_SAMPLE);
%         if fUseRaw
% %             mrWav_clu1 = single(median(tnWav_raw(:,:,viSpk_clu1), 3));
%             mrWav_clu1 = single(reshape(tnWav_raw(:,:,viSpk_clu1), nSamples, []));
%         else
% %             mrWav_clu1 = single(median(tnWav_spk(:,:,viSpk_clu1), 3));
%             mrWav_clu1 = single(reshape(tnWav_spk(:,:,viSpk_clu1), nSamples, []));
%         end
%     end
%     if fUseSd
%         mrPv1_clu(:,iClu) = std(mrWav_clu1,1,2);
%         if nPc>=2, mrPv2_clu(:,iClu) = mr_std2_(mrWav_clu1, nDelay)'; end
%         if nPc>=2, mrPv3_clu(:,iClu) = mr_std2_(mrWav_clu1, nDelay*2)'; end
%     else
%         [~, mrPv_clu1] = pca(mrWav_clu1, 'NumComponents', nPc);
%         mrPv1_clu(:,iClu) = mrPv_clu1(:,1);
%         if nPc>=2, mrPv2_clu(:,iClu) = mrPv_clu1(:,2); end
%         if nPc>=3, mrPv3_clu(:,iClu) = mrPv_clu1(:,3); end
%     end
% end
%
% if fUsePvCorr
%     func1 = @(x)1 - max(abs(xcorr_mr_(x, nShift)),[],3);
%     switch nPc
%         case 1
%             mrPcDist_clu = func1(mrPv1_clu);
%         case 2
%             mrPcDist_clu = (func1(mrPv1_clu) + func1(mrPv2_clu))/2;
%         case 3
%             mrPcDist_clu = (func1(mrPv1_clu) + func1(mrPv2_clu) + func1(mrPv3_clu))/3;
%     end
% else
%     [vrPc1_clu, vrPv1_clu] = pca(mrPv1_clu, 'NumComponents', 1);
%     switch nPc
%         case 1
%             mrPc_clu = vrPc1_clu;
%         case 2
%             [vrPc2_clu, vrPv2_clu] = pca(mrPv2_clu, 'NumComponents', 1);
%             mrPc_clu = [vrPc1_clu, vrPc2_clu];
%         case 3
%             [vrPc2_clu, vrPv2_clu] = pca(mrPv2_clu, 'NumComponents', 1);
%             [vrPc3_clu, vrPv3_clu] = pca(mrPv3_clu, 'NumComponents', 1);
%             mrPc_clu = [vrPc1_clu, vrPc2_clu, vrPc3_clu];
%     end
%     mrPcDist_clu = pdist2_(abs(mrPc_clu));
% end
% mrRealDist_clu = pdist2_(P.mrSiteXY(S_clu.viSite_clu,:));
%
% %mrPcDist_clu(sub2ind([nClu,nClu], 1:nClu, 1:nClu)) = nan;
% mrPcDist_clu(tril(true(nClu)) | mrRealDist_clu > MAX_REAL_DIST) = nan; %ignore bottom half
%
% % lower triangle only
% [vrMinDist_clu, viMinDist_clu] = min(mrPcDist_clu);
% vrRealDist_clu = mrRealDist_clu(sub2ind([nClu,nClu], 1:nClu, viMinDist_clu));
%
% vrMinDist_logz_clu = zeros(size(vrMinDist_clu));
% vi_NotNan = find(vrMinDist_clu > 0);
% vrMinDist_logz_clu(vi_NotNan) = zscore_(log(vrMinDist_clu(vi_NotNan)));
% % vrMinDist_logz_clu(2:end) = madscore_(log(vrMinDist_clu(2:end)));
%
% if nargout==0
%     figure; plot(vrMinDist_logz_clu, vrRealDist_clu, '.');
%     xlabel('min clu dist (log-MAD pc)'); ylabel('real clu dist (um)'); grid on;
%     vi_clu1 = find(vrMinDist_logz_clu < MAD_THRESH);
%     vi_clu2 = viMinDist_clu(vi_clu1);
%     vr_dist12 = vrMinDist_logz_clu(vi_clu1);
%     arrayfun(@(a,b,c)fprintf('(%d,%d,%0.2f), ', a,b,c), vi_clu1, vi_clu2, vr_dist12);
%     fprintf('\n');
% end
% end %func


%--------------------------------------------------------------------------
function [sd2, viRange1, viRange2] = mr_std2_(mr, nDelay)
    if nDelay==0, sd2 = std(mr,1,2); return; end

    % determine shift
    nT = size(mr,1);
    iShift1 = -round(nDelay/2);
    iShift2 = nDelay + iShift1;
    viRange1 = max((1:nT) + iShift1, 1);
    viRange2 = min((1:nT) + iShift2, nT);

    % viRange1 = (1:nT) + iShift1;
    % viRange2 = (1:nT) + iShift2;
    % vl12 = (viRange1>=1 & viRange1<=nT) & (viRange2>=1 & viRange2<=nT);
    % viRange1 = viRange1(vl12);
    % viRange2 = viRange2(vl12);

    mr1 = mr(viRange1,:);
    mr2 = mr(viRange2,:);

    sd2 = sqrt(abs(mean(mr1.*mr2,2) - mean(mr1,2).*mean(mr2,2)));
end %func
