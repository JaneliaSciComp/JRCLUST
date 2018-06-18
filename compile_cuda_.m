%--------------------------------------------------------------------------
% Compile CUDA codes for JRCLUST
% 10/5/17 JJJ: Error messages converted to warning
% 7/26/17 JJJ: Code cleanup and testing
function fSuccess = compile_cuda_(csFiles_cu)
    if nargin<1 || isempty(csFiles_cu)
        S_cfg = read_cfg_();
        csFiles_cu = S_cfg.csFiles_cu3; %version 3 cuda
    elseif ischar(csFiles_cu)
        csFiles_cu = {csFiles_cu};
    end
    t1 = tic;
    disp('Compiling CUDA codes...');
    fSuccess = 1;
    S_gpu = gpuDevice(1);
    if ispc()
        vcPath_nvcc = sprintf('"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v%0.1f\\bin\\nvcc"', S_gpu.ToolkitVersion);
    else
        vcPath_nvcc = '/usr/local/cuda/bin/nvcc';
    end

    % C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.5\bin
    for i=1:numel(csFiles_cu)
        vcFile_ = jrcpath_(csFiles_cu{i});
        vcCmd1 = sprintf('%s -ptx -m 64 -arch sm_35 "%s"', vcPath_nvcc, vcFile_);
        fprintf('\t%s\n\t', vcCmd1);
        try
            status = system(vcCmd1);
            fSuccess = fSuccess && (status==0);
        catch
            fprintf('\tWarning: CUDA could not be compiled: %s\n', vcFile_);
        end
    end
    if ~fSuccess
        fprintf('\tWarning: CUDA could not be compiled but JRCLUST may work fine. If not, install CUDA toolkit v%0.1f and run "jrc install".\n', S_gpu.ToolkitVersion);
    end
    fprintf('\tFinished compiling, took %0.1fs\n', toc(t1));
end %func
