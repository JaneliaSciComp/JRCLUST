%--------------------------------------------------------------------------
% 2017/12/5 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_2_(S_clu, P, viClu_update) % works for only the latest format
    % use extra padding for spkRaw to time-shift nShift

    P.useGPU = 0;

    if nargin<3, viClu_update = []; end
    if ~isfield(S_clu, 'mrWavCor'), viClu_update = []; end
    fprintf('Computing waveform correlation...'); t1 = tic;
    ctmrWav_clu = {S_clu.tmrWav_raw_clu, S_clu.tmrWav_raw_lo_clu, S_clu.tmrWav_raw_hi_clu};

    nClu = S_clu.nClusters;
    mrWavCor = gpuArray_(zeros(nClu), P.useGPU);
    if isempty(viClu_update)
        vlClu_update = true(nClu, 1);
        mrWavCor0 = [];
    else
        vlClu_update = false(nClu, 1);
        vlClu_update(viClu_update) = 1;
        mrWavCor0 = S_clu.mrWavCor;
        nClu_pre = size(mrWavCor0, 1);
        vlClu_update((1:nClu) > nClu_pre) = 1;
    end
    [cviShift1, cviShift2] = calc_shift_range_(P);
    cell_5args = {P, vlClu_update, mrWavCor0, cviShift1, cviShift2};
    try
        parfor iClu2 = 1:nClu %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, S_clu.viSite_clu, iClu2, cell_5args);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    catch
        fprintf('S_clu_wavcor_: parfor failed. retrying for loop\n');
        for iClu2 = 1:nClu %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, S_clu.viSite_clu, iClu2, cell_5args);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    end

    mrWavCor = max(mrWavCor, mrWavCor'); %make it symmetric
    mrWavCor(mrWavCor==0) = nan;
    mrWavCor = gather_(mrWavCor);
    fprintf('\ttook %0.1fs\n', toc(t1));
end %func
