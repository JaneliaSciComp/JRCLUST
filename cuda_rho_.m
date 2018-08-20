%--------------------------------------------------------------------------
function vrRho1 = cuda_rho_(mrFet12, viiSpk12_ord, n1, n2, dc2, P)
    % Ultimately use CUDA to do this distance computation
    % mrFet12_: already in GPU
    % viiSpk12: ordered list
    persistent CK nC_
    if nargin==0, nC_ = 0; return; end
    if isempty(nC_), nC_ = 0; end
    [nC, n12] = size(mrFet12); %nc is constant with the loop
    % if getOr(P, 'f_dpclus', 1)
    dn_max = int32(round((n1+n2) / P.nTime_clu));
    % else
    %     dn_max = int32(round((n1+n2)));
    % end
    nC_max = getOr(P, 'nC_max', 45);
    dc2 = single(dc2);
    FLAG_FIXN = 0; %flag variable number of neighbors (otherwise fixed to 2*dn_max+1)
    if P.useGPU && FLAG_FIXN == 0
        try
            if (nC_ ~= nC) % create cuda kernel
                nC_ = nC;
                CK = parallel.gpu.CUDAKernel('jrc_cuda_rho.ptx','jrc_cuda_rho.cu');
                CK.ThreadBlockSize = [P.nThreads, 1];
                CK.SharedMemorySize = 4 * P.CHUNK * (2 + nC_max + 2 * P.nThreads); % @TODO: update the size
            end
            CK.GridSize = [ceil(n1 / P.CHUNK / P.CHUNK), P.CHUNK]; %MaxGridSize: [2.1475e+09 65535 65535]
            vrRho1 = zeros([1, n1], 'single', 'gpuArray');
            vnConst = int32([n1, n12, nC, dn_max, getOr(P, 'fDc_spk', 0)]);
            vrRho1 = feval(CK, vrRho1, mrFet12, viiSpk12_ord, vnConst, dc2);
            return;
        catch
            disperr_('CUDA kernel failed. Re-trying in CPU.');
            nC_ = 0;
            %         fprintf(2, 'CUDA kernel failed. Re-trying in CPU\n');
        end
    end

    viiSpk1_ord= viiSpk12_ord(1:n1)';
    mlKeep12_ = abs(bsxfun(@minus, viiSpk12_ord, viiSpk1_ord)) <= dn_max;
    if FLAG_FIXN == 0
        vrRho1 = sum(eucl2_dist_(mrFet12, mrFet12(:,1:n1)) < dc2 & mlKeep12_) - 1; %do not include self
        vrRho1 = single(vrRho1 ./ sum(mlKeep12_));
    else
        mlKeep12_(1:min(2*dn_max+1,n12), viiSpk1_ord <= min(dn_max,n1)) = 1;
        mlKeep12_(max(n12-2*dn_max,1):n12,  viiSpk1_ord >= max(n1-dn_max,1)) = 1;
        vrRho1 = sum(eucl2_dist_(mrFet12, mrFet12(:,1:n1)) < dc2 & mlKeep12_) - 1; %do not include self
        vrRho1 = single(vrRho1 ./ (2*dn_max));
    end
end
