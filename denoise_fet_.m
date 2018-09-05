%--------------------------------------------------------------------------
% 10/12/17 JJJ: site by site denoising
function spikeFeatures_ = denoise_fet_(spikeFeatures, P, vlRedo_spk)
    % denoise_fet_() to reset
    % cluster based averaging fet cleanup
    % set repeat parameter?
    nRepeat_fet = getOr(P, 'nRepeat_fet', 0);
    if nRepeat_fet==0, spikeFeatures_ = spikeFeatures; return ;end

    S0 = get(0, 'UserData');
    fprintf('Denoising features using nneigh\n\t'); t1=tic;
    nSites = numel(P.chanMap);
    nC = size(spikeFeatures,1);
    try
        nC_max = getOr(P, 'nC_max', 45);
        CK = parallel.gpu.CUDAKernel('jrc_cuda_nneigh.ptx','jrc_cuda_nneigh.cu');
        CK.ThreadBlockSize = [P.nThreads, 1];
        CK.SharedMemorySize = 4 * P.CHUNK * (1 + nC_max + 2 * P.nThreads); % @TODO: update the size
    catch
        disperr_('denoise_fet_: Cuda init error');
    end
    try
        for iSite = 1:nSites
            viSpk1 = S0.cviSpk_site{iSite};
            if isempty(viSpk1), continue; end
            if ~isempty(vlRedo_spk)
                viSpk1 = viSpk1(vlRedo_spk(viSpk1));
                if isempty(viSpk1), continue; end
            end
            trFet12 = permute(spikeFeatures(:,:,viSpk1), [1,3,2]);
            mrFet1 = gpuArray(trFet12(:,:,1));
            mrFet2 = gpuArray(trFet12(:,:,2));
            viiSpk1_ord = gpuArray(rankOrder(viSpk1, 'ascend'));
            n1 = numel(viSpk1);
            dn_max = int32(round(n1 / P.nTime_clu));
            vrDelta1 = zeros([1, n1], 'single', 'gpuArray');
            viNneigh1 = zeros([1, n1], 'uint32', 'gpuArray');
            vnConst = int32([n1, nC, dn_max]);
            CK.GridSize = [ceil(n1 / P.CHUNK / P.CHUNK), P.CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]

            for iRepeat = 1:nRepeat_fet
                [vrDelta1, viNneigh1] = feval(CK, vrDelta1, viNneigh1, mrFet1, viiSpk1_ord, vnConst);
                mrFet1 = (mrFet1 + mrFet1(:,viNneigh1)) / 2;
                mrFet2 = (mrFet2 + mrFet2(:,viNneigh1)) / 2;
                %             mrD11 = set_diag_(eucl2_dist_(mrFet1, mrFet1), nan);
                %               [vrDelta2, viNneigh2] = min(mrD11);
            end

            spikeFeatures_(:,1,viSpk1) = gather_(mrFet1);
            spikeFeatures_(:,2,viSpk1) = gather_(mrFet2);
            [mrFet1, mrFet2, vrDelta1, viNneigh1] = deal([]);
            fprintf('.');
        end
    catch
        disperr_('denoise_fet_: CUDA init error');
    end
    % [vrRho, vrDc2_site] = gather_(vrRho, vrDc2_site);
    fprintf('\n\ttook %0.1fs\n', toc(t1));
end
