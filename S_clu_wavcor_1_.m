%--------------------------------------------------------------------------
% 6/29/17 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_1_(S_clu, P, viClu_update)
    % symmetric matrix and common basis comparison only

    fUsePeak2 = 0;
    nInterp_merge = get_set_(P, 'nInterp_merge', 1); % set to 1 to disable
    fDrift_merge = get_set_(P, 'fDrift_merge', 0);
    P.useGPU = 0;
    nShift = ceil(P.spkRefrac_ms / 1000 * P.sRateHz * nInterp_merge); % +/-n number of samples to compare time shift
    % nShift = 0;
    fWaveform_raw = 0; % get_set_(P, 'fWavRaw_merge', 1); TW

    fZeroStart_raw = get_set_(P, 'fZeroStart_raw', 0);
    fRankCorr_merge = get_set_(P, 'fRankCorr_merge', 0);
    fMode_cor = 1; %0: pearson, 1:no mean subt pearson

    if nargin<3, viClu_update = []; end
    if ~isfield(S_clu, 'mrWavCor'), viClu_update = []; end
    fprintf('Computing waveform correlation...'); t1 = tic;
    if fWaveform_raw
        tmrWav_clu = trim_spkraw_(S_clu.tmrWav_raw_clu, P);
    else
        tmrWav_clu = S_clu.tmrWav_spk_clu;
    end
    % tmrWav_clu = ifeq_(fWaveform_raw, S_clu.tmrWav_raw_clu, S_clu.tmrWav_spk_clu);
    % tmrWav_clu = gpuArray_(tmrWav_clu, P.useGPU);
    if fDrift_merge && fWaveform_raw % only works on raw
        tmrWav_lo_clu = trim_spkraw_(S_clu.tmrWav_raw_lo_clu, P);
        tmrWav_hi_clu = trim_spkraw_(S_clu.tmrWav_raw_hi_clu, P);
        ctmrWav_clu = {tmrWav_clu, tmrWav_lo_clu, tmrWav_hi_clu};
    else
        ctmrWav_clu = {tmrWav_clu};
    end
    if fZeroStart_raw
        ctmrWav_clu = cellfun(@(x)zero_start_(x), ctmrWav_clu, 'UniformOutput', 0);
    end
    if fRankCorr_merge
        ctmrWav_clu = cellfun(@(x)rankorder_mr_(x,0), ctmrWav_clu, 'UniformOutput', 0);
    end
    if nInterp_merge>1
        ctmrWav_clu = cellfun(@(x)interpft_(x, nInterp_merge), ctmrWav_clu, 'UniformOutput', 0);
    end
    nClu = S_clu.nClu;
    if fUsePeak2
        [viSite_clu, viSite2_clu, viSite3_clu] = S_clu_peak2_(S_clu);
        cviSite_clu = {viSite_clu, viSite2_clu, viSite3_clu};
    else
        viSite_clu = S_clu.viSite_clu;
        cviSite_clu = {viSite_clu};
    end
    mrWavCor = gpuArray_(zeros(nClu), P.useGPU);
    nSites_spk = P.maxSite*2+1-P.nSites_ref;
    nT = size(tmrWav_clu, 1);
    [cviShift1, cviShift2] = shift_range_(gpuArray_(int32(nT), P.useGPU), nShift);
    % viLags = 1:numel(cviShift1);
    if isempty(viClu_update)
        vlClu_update = true(nClu, 1);
        mrWavCor0 = [];
        fParfor = 1;
    else
        fParfor = 0;
        vlClu_update = false(nClu, 1);
        vlClu_update(viClu_update) = 1;
        mrWavCor0 = S_clu.mrWavCor;
        nClu_pre = size(mrWavCor0, 1);
        vlClu_update((1:nClu) > nClu_pre) = 1;
    end
    cell_5args = {vlClu_update, cviShift1, cviShift2, mrWavCor0, fMode_cor};
    if fParfor
        try
            parfor iClu2 = 1:nClu %parfor speedup: 4x
                vrWavCor2 = clu_wavcor_(ctmrWav_clu, cviSite_clu, P, cell_5args, iClu2);
                if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
            end
        catch
            fprintf('S_clu_wavcor_: parfor failed. retrying for loop\n');
            fParfor = 0;
        end
    end
    if ~fParfor
        for iClu2 = 1:nClu %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_(ctmrWav_clu, cviSite_clu, P, cell_5args, iClu2);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    end

    mrWavCor = max(mrWavCor, mrWavCor'); %make it symmetric
    mrWavCor(mrWavCor==0) = nan;
    mrWavCor = gather_(mrWavCor);
    fprintf('\ttook %0.1fs\n', toc(t1));
end %func
