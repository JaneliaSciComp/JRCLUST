%--------------------------------------------------------------------------
function [vrDelta1, viNneigh1] = cuda_delta_(mrFet12, viiSpk12_ord, viiRho12_ord, n1, n2, dc2, P)
    % Ultimately use CUDA to do this distance computation
    persistent CK nC_
    if nargin==0, nC_ = 0; return; end
    if isempty(nC_), nC_ = 0; end
    [nC, n12] = size(mrFet12); %nc is constant with the loop
    dn_max = int32(round((n1+n2) / P.nTime_clu));
    nC_max = get_set_(P, 'nC_max', 45);
    if P.useGPU
        try
            if (nC_ ~= nC) % create cuda kernel
                nC_ = nC;
                CK = parallel.gpu.CUDAKernel('jrc_cuda_delta.ptx','jrc_cuda_delta.cu');
                CK.ThreadBlockSize = [P.nThreads, 1];
                CK.SharedMemorySize = 4 * P.CHUNK * (3 + nC_max + 2*P.nThreads); % @TODO: update the size
            end
            CK.GridSize = [ceil(n1 / P.CHUNK / P.CHUNK), P.CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]
            vrDelta1 = zeros([1, n1], 'single', 'gpuArray');
            viNneigh1 = zeros([1, n1], 'uint32', 'gpuArray');
            vnConst = int32([n1, n12, nC, dn_max, get_set_(P, 'fDc_spk', 0)]);
            [vrDelta1, viNneigh1] = feval(CK, vrDelta1, viNneigh1, mrFet12, viiSpk12_ord, viiRho12_ord, vnConst, dc2);
            % [vrDelta1_, viNneigh1_] = deal(vrDelta1, viNneigh1);
            return;
        catch
            disperr_('CUDA kernel failed. Re-trying in CPU.');
            nC_ = 0;
        end
    end

    mrDist12_ = eucl2_dist_(mrFet12, mrFet12(:,1:n1)); %not sqrt
    mlRemove12_ = bsxfun(@ge, viiRho12_ord, viiRho12_ord(1:n1)') ...
    | abs(bsxfun(@minus, viiSpk12_ord, viiSpk12_ord(1:n1)')) > dn_max;
    mrDist12_(mlRemove12_) = nan;
    [vrDelta1, viNneigh1] = min(mrDist12_);
end
