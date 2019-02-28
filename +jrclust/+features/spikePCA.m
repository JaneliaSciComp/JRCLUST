function [features1, features2, features3] = spikePCA(spikeWindows, hCfg)
    %SPIKEPCA Project spikes onto their principal components
    % if strcmpi(hCfg.vcFet, 'pca')
    prVecs = spikePrVecs(spikeWindows, hCfg);
    % else
    %     mrPv_global = get0_('mrPv_global');
    %     if isempty(mrPv_global)
    %         [mrPv_global, vrD_global] = spikePrVecs(spikeWindows, hCfg);
    %         [mrPv_global, vrD_global] = jrclust.utils.tryGather(mrPv_global, vrD_global);
    %         set0_(mrPv_global, vrD_global);
    %     end

    %     prVecs = mrPv_global;
    % end

    [features1, features2, features3] = projectInterp(spikeWindows, prVecs, hCfg);
end

%% LOCAL FUNCTIONS
function [prVecs, eigVals] = spikePrVecs(spikeWindows, hCfg)
    %SPIKEPRVECS Get principal vectors for spikes
    MAX_SAMPLE = 10000;

    % subsample spikes on their primary sites up to MAX_SAMPLE
    ss = jrclust.utils.subsample(1:size(spikeWindows, 2), MAX_SAMPLE);
    spikeSample = spikeWindows(:, ss, 1);

    covMat = (spikeSample*spikeSample')/(numel(ss)-1);

    [eigVecs, eigVals] = eig(covMat);
    % MATLAB returns value-vector pairs smallest to largest; flip them
    eigVecs = jrclust.utils.zscore(fliplr(eigVecs));
    eigVals = flipud(diag(eigVals));

    % spike center should be negative
    iMid = 1-hCfg.evtWindowSamp(1); % sample where the spike is said to occur
    sgn = (eigVecs(iMid, :) < 0) * 2 - 1; % 1 or -1 depending on the sign
    prVecs = bsxfun(@times, eigVecs, sgn);
end

function [features1, features2, features3] = projectInterp(spikeWindows, prVecs, hCfg)
    %PROJECTINTERP Project spikes onto principal vectors
    shape = size(spikeWindows);
    if ismatrix(spikeWindows)
        shape(end + 1) = 1;
    end

    % lay the deck out left to right
    spikeWaveforms = reshape(spikeWindows, shape(1), []);

    prVecs = jrclust.utils.tryGather(prVecs);
    % project spikes onto first PC
    features1 = reshape(prVecs(:, 1)'*spikeWaveforms, shape(2:3))';

    if hCfg.nPCsPerSite >= 2 % project spikes onto second PC
        features2 = reshape(prVecs(:, 2)'*spikeWaveforms, shape(2:3))';
    else
        features2 = [];
    end
    if hCfg.nPCsPerSite == 3 % project spikes onto third PC
        features3 = reshape(prVecs(:, 3)'*spikeWaveforms, shape(2:3))';
    else
        features3 = [];
    end

    if ~hCfg.interpPC
        return;
    end

    % find optimal delay by interpolating
    vr1 = prVecs(:, 1);
    vi0 = (1:numel(vr1))';
    shifts = [0, -1, -.5, .5, 1];

    mrPv1 = zeros(numel(vr1), numel(shifts), 'like', prVecs);
    mrPv1(:, 1) = vr1;

    for iShift = 2:numel(shifts)
        mrPv1(:, iShift) = zscore(interp1(vi0, vr1, vi0+shifts(iShift), 'pchip', 'extrap'));
    end

    % find shift that maximizes the projection
    [~, viMax_spk] = max(abs(mrPv1'*spikeWindows(:, :, 1)));

    for iShift=2:numel(shifts)
        viSpk2 = find(viMax_spk == iShift);
        if isempty(viSpk2)
            continue;
        end

        mrWav_spk2 = reshape(spikeWindows(:, viSpk2, :), shape(1), []);
        features1(:, viSpk2) = reshape(mrPv1(:, iShift)'*mrWav_spk2, [], shape(3))';
    end
end
