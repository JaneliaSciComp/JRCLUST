function success = compileCUDA(nvccPath)
    %COMPILECUDA Compile CUDA codes for JRCLUST
    gpuD = gpuDevice(1);
    if nargin < 1
        if ispc()
            nvccPath = sprintf('"C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v%0.1f\\bin\\nvcc.exe"', gpuD.ToolkitVersion);
        else
            [ecode, out] = system('which nvcc');
            if ecode == 0 % successfully found nvcc on path
                nvccPath = strip(out);
            else % fall back to default
                nvccPath = '/usr/local/cuda/bin/nvcc';
            end
        end
    end
    basedir = fullfile(jrclust.utils.basedir(), '+jrclust', '+CUDA');
    cudaFiles = dir(fullfile(basedir, '*.cu'));
    cudaFiles = {cudaFiles.name};

    t1 = tic;
    disp('Compiling CUDA codes...');
    success = 1;

    for i = 1:numel(cudaFiles)
        iFileCU = fullfile(basedir, cudaFiles{i});
        iFilePTX = fullfile(basedir, strrep(cudaFiles{i}, '.cu', '.ptx'));
        cmd = sprintf('%s -ptx -m 64 -arch sm_35 "%s" --output-file "%s"', nvccPath, iFileCU, iFilePTX);

        fprintf('\t%s\n\t', cmd);

        try
            status = system(cmd);
            success = success && (status == 0);
        catch ME
            warning('Could not compile %s: %s\n', iFileCU, ME.message);
        end
    end

    if ~success
        warning('CUDA could not be compiled but JRCLUST may work fine. If not, install CUDA toolkit v%0.1f and run "jrc compile".\n', gpuD.ToolkitVersion);
    end

    fprintf('Finished compiling, took %0.1fs\n', toc(t1));
end
