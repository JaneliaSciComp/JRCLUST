function mrWavCor = S_clu_wavcor_(hClust, hCfg, updateMe)
    if nargin < 3
        updateMe = [];
    end

    fUse_old = true; % TODO: investigate S_clu_wavcor_2_
    if hCfg.spkLim_factor_merge < hCfg.evtWindowRawFactor && ~fUse_old
        mrWavCor = S_clu_wavcor_2_(hClust, hCfg, updateMe); % newer format v3.1.8
    else
        mrWavCor = S_clu_wavcor_1_(hClust, hCfg, updateMe); % older format
    end
end

%% LOCAL FUNCTIONS
function mrWavCor = S_clu_wavcor_1_(hClust, hCfg, updateMe)
    % symmetric matrix and common basis comparison only
    if nargin < 3 || isempty(hClust.mrWavCor)
        updateMe = [];
    end

    hCfg.useGPU = false;

    % +/-n number of samples to compare time shift
    nShift = ceil(hCfg.sampleRate*hCfg.nInterp_merge*hCfg.refracIntms/1000);

    % not in default param set, but can be overridden if you really want
    fUsePeak2 = hCfg.getOr('fUsePeak2', false);
    fWaveform_raw = hCfg.getOr('fWavRaw_merge', true); % revert TW: get_set_(P, 'fWavRaw_merge', 0);
    fZeroStart_raw = hCfg.getOr('fZeroStart_raw', false);
    fRankCorr_merge = hCfg.getOr('fRankCorr_merge', false);
    fMode_cor = hCfg.getOr('fMode_cor', 1); % 0: pearson, 1: no mean subt pearson

    if hCfg.verbose
        fprintf('Computing waveform correlation...');
        t1 = tic;
    end

    if fWaveform_raw
        meanWfGlobal = trim_spkraw_(hClust.meanWfGlobalRaw, hCfg);
    else
        meanWfGlobal = hClust.meanWfGlobal;
    end

    if hCfg.fDrift_merge && fWaveform_raw % only works on raw
        meanWfRawLow = trim_spkraw_(hClust.meanWfRawLow, hCfg);
        meanWfRawHigh = trim_spkraw_(hClust.meanWfRawHigh, hCfg);
        meanWfSet = {meanWfGlobal, meanWfRawLow, meanWfRawHigh};
    else
        meanWfSet = {meanWfGlobal};
    end

    if fZeroStart_raw
        meanWfSet = cellfun(@(x) zeroStart(x), meanWfSet, 'UniformOutput', 0);
    end

    if fRankCorr_merge
        meanWfSet = cellfun(@(x) rankorder_mr_(x), meanWfSet, 'UniformOutput', 0);
    end

    if hCfg.nInterp_merge > 1
        meanWfSet = cellfun(@(x) interpft_(x, hCfg.nInterp_merge), meanWfSet, 'UniformOutput', 0);
    end

    if fUsePeak2
        [peaks1, peaks2, peaks3] = getSecondaryPeaks(hClust);
        clusterSites = {peaks1, peaks2, peaks3};
    else
        peaks1 = hClust.clusterSites;
        clusterSites = {peaks1};
    end

    mrWavCor = jrclust.utils.tryGpuArray(zeros(hClust.nClusters), hCfg.useGPU);

    nT = size(meanWfGlobal, 1);
    [cviShift1, cviShift2] = jrclust.utils.shiftRange(jrclust.utils.tryGpuArray(int32(nT), hCfg.useGPU), nShift);

    if isempty(updateMe)
        vlClu_update = true(hClust.nClusters, 1);
        mrWavCor0 = [];
        fParfor = 1;
    else
        fParfor = 0;
        vlClu_update = false(hClust.nClusters, 1);
        vlClu_update(updateMe) = 1;
        mrWavCor0 = hClust.mrWavCor;
        nClu_pre = size(mrWavCor0, 1);
        vlClu_update((1:hClust.nClusters) > nClu_pre) = 1;
    end
    cell_5args = {vlClu_update, cviShift1, cviShift2, mrWavCor0, fMode_cor};

    if fParfor
        try
            parfor iClu2 = 1:hClust.nClusters %parfor speedup: 4x
                vrWavCor2 = clu_wavcor_(meanWfSet, clusterSites, hCfg, cell_5args, iClu2);
                if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
            end
        catch
            fprintf('S_clu_wavcor_: parfor failed. retrying for loop\n');
            fParfor = 0;
        end
    end

    if ~fParfor
        for iClu2 = 1:hClust.nClusters %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_(meanWfSet, clusterSites, hCfg, cell_5args, iClu2);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    end

    mrWavCor = max(mrWavCor, mrWavCor'); %make it symmetric
    mrWavCor(mrWavCor==0) = nan;
    mrWavCor = jrclust.utils.tryGather(mrWavCor);

    if hCfg.verbose
        fprintf('\ttook %0.1fs\n', toc(t1));
    end
end %func

%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function vrWavCor2 = clu_wavcor_(ctmrWav_clu, cviSite_clu, P, cell_5args, iClu2)

    [vlClu_update, cviShift1, cviShift2, mrWavCor0, fMode_cor] = deal(cell_5args{:});
    if numel(cviSite_clu) == 1
        viSite_clu = cviSite_clu{1};
        fUsePeak2 = 0;
    else
        [viSite_clu, viSite2_clu, viSite3_clu] = deal(cviSite_clu{:});
        fUsePeak2 = 1;
    end
    nClusters = numel(viSite_clu);
    iSite_clu2 = viSite_clu(iClu2);
    if iSite_clu2==0 || isnan(iSite_clu2), vrWavCor2 = []; return; end
    viSite2 = P.miSites(:,iSite_clu2);
    % if fMaxSite_excl, viSite2 = viSite2(2:end); end
    if fUsePeak2
        viClu1 = find(viSite_clu == iSite_clu2 | viSite2_clu == iSite_clu2 | viSite3_clu == iSite_clu2 | ...
        viSite_clu == viSite2_clu(iClu2) | viSite_clu == viSite3_clu(iClu2)); %viSite2_clu == viSite2_clu(iClu2)
    else
        %     maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
        maxDist_site_um = get_set_(P, 'maxDist_site_merge_um', 35);
        viClu1 = find(ismember(viSite_clu, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));
    end

    vrWavCor2 = zeros(nClusters, 1, 'single');
    viClu1(viClu1 >= iClu2) = []; % symmetric matrix comparison
    if isempty(viClu1), return; end
    cmrWav_clu2 = cellfun(@(x)x(:,viSite2,iClu2), ctmrWav_clu, 'UniformOutput', 0);
    ctmrWav_clu1 = cellfun(@(x)x(:,viSite2,viClu1), ctmrWav_clu, 'UniformOutput', 0);
    for iClu11 = 1:numel(viClu1)
        iClu1 = viClu1(iClu11);
        if ~vlClu_update(iClu1) && ~vlClu_update(iClu2)
            vrWavCor2(iClu1) = mrWavCor0(iClu1, iClu2);
        else
            iSite_clu1 = viSite_clu(iClu1);
            if iSite_clu1==0 || isnan(iSite_clu1), continue; end
            if iSite_clu1 == iSite_clu2
                cmrWav_clu2_ = cmrWav_clu2;
                cmrWav_clu1_ = cellfun(@(x)x(:,:,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            else
                viSite1 = P.miSites(:,iSite_clu1);
                viSite12 = find(ismember(viSite2, viSite1));
                if isempty(viSite12), continue; end
                cmrWav_clu2_ = cellfun(@(x)x(:,viSite12), cmrWav_clu2, 'UniformOutput', 0);
                cmrWav_clu1_ = cellfun(@(x)x(:,viSite12,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            end
            vrWavCor2(iClu1) = maxCor_drift_(cmrWav_clu2_, cmrWav_clu1_, cviShift1, cviShift2, fMode_cor);
        end
    end %iClu2 loop
end %func

%--------------------------------------------------------------------------
% 10/23/17 JJJ: find max correlation pair (combining drift and temporal shift)
function maxCor = maxCor_drift_(cmr1, cmr2, cviShift1, cviShift2, fMode_cor)
    if nargin < 5
        fMode_cor = 0;
    end % pearson corr

    assert_(numel(cmr1) == numel(cmr2), 'maxCor_drift_: numel must be the same');

    P = get0_('P');

    if numel(cmr1) == 1
        if strcmpi(get_set_(P, 'autoMergeCriterion', 'xcorr'), 'dist')
            a_=reshape(cmr1{1},[],1); % TW?
            b_=reshape(cmr2{1},[],1); % TW?

            % m_=max([var(a_) var(b_)]); % TW?
            % tmp_=cov(a_,b_)/m_; % TW?
            % c_(counter) = tmp_(1,2); % TW?
            maxCor= 1-sum(abs(a_-b_))/max([sum(abs(a_));sum(abs(b_))]); % TW?
        else % default behavior: Pearson coeff
            maxCor = max(xcorr2_mr_(cmr1{1}, cmr2{1}, cviShift1, cviShift2));
        end
    else
        if strcmpi(get_set_(P, 'autoMergeCriterion', 'xcorr'), 'dist')
            % begin TW block
            c_ = zeros(size(tr1, 3), 1);
            for counter=1:size(tr1,3)
                a_=reshape(tr1(:,:,counter),[],1);
                b_=reshape(tr2(:,:,counter),[],1);

                % m_=max([var(a_) var(b_)]);
                % tmp_=cov(a_,b_)/m_;
                % c_(counter) = tmp_(1,2);
                c_(counter) = -sum(abs(a-b))/max(abs([a_ b_]));
            end

            maxCor=max(c_);
            % end TW block
        else
            tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
            tr2 = cat(3, cmr2{:});
            nDrift = numel(cmr1);
            nShift = numel(cviShift1);
            vrCor = zeros(1, nShift);
            for iShift = 1:nShift
                mr1_ = reshape(tr1(cviShift1{iShift},:,:), [], nDrift);
                mr2_ = reshape(tr2(cviShift2{iShift},:,:), [], nDrift);
                if fMode_cor == 0
                    mr1_ = bsxfun(@minus, mr1_, sum(mr1_)/size(mr1_,1));
                    mr2_ = bsxfun(@minus, mr2_, sum(mr2_)/size(mr2_,1));
                    mr1_ = bsxfun(@rdivide, mr1_, sqrt(sum(mr1_.^2)));
                    mr2_ = bsxfun(@rdivide, mr2_, sqrt(sum(mr2_.^2)));
                    vrCor(iShift) = max(max(mr1_' * mr2_));
                else
                    mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
                    vrCor(iShift) = max(mr12_(:));
                end % if
            end % for

            maxCor = max(vrCor);
        end % if

        %     tic
        %     mrCor = zeros(numel(cmr1), numel(cmr2));
        %     for i1=1:numel(cmr1)
        %         for i2=1:numel(cmr2)
        %             mrCor(i1,i2) = max(xcorr2_mr_(cmr1{i1}, cmr2{i2}, cviShift1, cviShift2));
        %         end
        %     end
        %     maxCor = max(mrCor(:));
        %     toc
    end % if
end

%--------------------------------------------------------------------------
% 2017/12/5 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_2_(S_clu, P, viClu_update) % works for only the latest format
    % use extra padding for spkRaw to time-shift nShift

    P.fGpu = 0;

    if nargin<3, viClu_update = []; end
    if ~isfield(S_clu, 'mrWavCor'), viClu_update = []; end
    fprintf('Computing waveform correlation...'); t1 = tic;
    ctmrWav_clu = {S_clu.tmrWav_raw_clu, S_clu.tmrWav_raw_lo_clu, S_clu.tmrWav_raw_hi_clu};

    nClusters = S_clu.nClusters;
    mrWavCor = jrclust.utils.tryGpuArray(zeros(nClusters), P.fGpu);
    if isempty(viClu_update)
        vlClu_update = true(nClusters, 1);
        mrWavCor0 = [];
    else
        vlClu_update = false(nClusters, 1);
        vlClu_update(viClu_update) = 1;
        mrWavCor0 = S_clu.mrWavCor;
        nClu_pre = size(mrWavCor0, 1);
        vlClu_update((1:nClusters) > nClu_pre) = 1;
    end
    [cviShift1, cviShift2] = calc_shift_range_(P);
    cell_5args = {P, vlClu_update, mrWavCor0, cviShift1, cviShift2};
    try
        parfor iClu2 = 1:nClusters %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, S_clu.viSite_clu, iClu2, cell_5args);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    catch
        fprintf('S_clu_wavcor_: parfor failed. retrying for loop\n');
        for iClu2 = 1:nClusters %parfor speedup: 4x
            vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, S_clu.viSite_clu, iClu2, cell_5args);
            if ~isempty(vrWavCor2), mrWavCor(:, iClu2) = vrWavCor2; end
        end
    end

    mrWavCor = max(mrWavCor, mrWavCor'); %make it symmetric
    mrWavCor(mrWavCor==0) = nan;
    mrWavCor = jrclust.utils.tryGather(mrWavCor);
    fprintf('\ttook %0.1fs\n', toc(t1));
end %func

%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, viSite_clu, iClu2, cell_5args)
    [P, vlClu_update, mrWavCor0, cviShift1, cviShift2] = deal(cell_5args{:});
    nClusters = numel(viSite_clu);
    iSite_clu2 = viSite_clu(iClu2);
    if iSite_clu2==0 || isnan(iSite_clu2), vrWavCor2 = []; return; end
    viSite2 = P.miSites(:,iSite_clu2);
    maxDist_site_um = get_set_(P, 'maxDist_site_merge_um', 35);
    viClu1 = find(ismember(viSite_clu, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));

    vrWavCor2 = zeros(nClusters, 1, 'single');
    viClu1(viClu1 >= iClu2) = []; % symmetric matrix comparison
    if isempty(viClu1), return; end
    cmrWav_clu2 = cellfun(@(x)x(:,viSite2,iClu2), ctmrWav_clu, 'UniformOutput', 0);
    ctmrWav_clu1 = cellfun(@(x)x(:,viSite2,viClu1), ctmrWav_clu, 'UniformOutput', 0);
    for iClu11 = 1:numel(viClu1)
        iClu1 = viClu1(iClu11);
        if ~vlClu_update(iClu1) && ~vlClu_update(iClu2)
            vrWavCor2(iClu1) = mrWavCor0(iClu1, iClu2);
        else
            iSite_clu1 = viSite_clu(iClu1);
            if iSite_clu1==0 || isnan(iSite_clu1), continue; end
            if iSite_clu1 == iSite_clu2
                cmrWav_clu2_ = cmrWav_clu2;
                cmrWav_clu1_ = cellfun(@(x)x(:,:,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            else
                viSite1 = P.miSites(:,iSite_clu1);
                viSite12 = find(ismember(viSite2, viSite1));
                if isempty(viSite12), continue; end
                cmrWav_clu2_ = cellfun(@(x)x(:,viSite12), cmrWav_clu2, 'UniformOutput', 0);
                cmrWav_clu1_ = cellfun(@(x)x(:,viSite12,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            end
            vrWavCor2(iClu1) = maxCor_drift_2_(cmrWav_clu2_, cmrWav_clu1_, cviShift1, cviShift2);
        end
    end %iClu2 loop
end %func

% 10/23/17 JJJ: find max correlation pair (combining drift and temporal shift)
%--------------------------------------------------------------------------
function maxCor = maxCor_drift_2_(cmr1, cmr2, cviShift1, cviShift2)
    % assert(numel(cmr1) == numel(cmr2), 'maxCor_drift_: numel must be the same');
    nDrift = numel(cmr1);
    if nDrift == 1
        maxCor = max(xcorr2_mr_(cmr1{1}, cmr2{1}, cviShift1, cviShift2));
    else
        tr1 = cat(3, cmr1{:}); %nT x nC x nDrifts
        tr2 = cat(3, cmr2{:});
        nShift = numel(cviShift1);
        vrCor = zeros(1, nShift);
        for iShift = 1:nShift
            mr1_ = reshape(jrclust.utils.meanSubtract(tr1(cviShift1{iShift},:,:)), [], nDrift);
            mr2_ = reshape(jrclust.utils.meanSubtract(tr2(cviShift2{iShift},:,:)), [], nDrift);

            mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
            vrCor(iShift) = max(mr12_(:));
        end
        maxCor = max(vrCor);
    end
end

function meanWf = trim_spkraw_(meanWf, hCfg)
    nSamples_raw = diff(hCfg.evtWindowRawSamp) + 1;
    spkLim_merge = round(hCfg.evtWindowSamp * hCfg.spkLim_factor_merge);
    nSamples_raw_merge = diff(spkLim_merge) + 1;

    if size(meanWf, 1) <= nSamples_raw_merge
        return;
    end

    lim_merge = [spkLim_merge(1) - hCfg.evtWindowRawSamp(1) + 1,  nSamples_raw - (hCfg.evtWindowRawSamp(2) - spkLim_merge(2))];
    meanWf = meanWf(lim_merge(1):lim_merge(2), :, :);
    meanWf = jrclust.utils.meanSubtract(meanWf);
end

function meanWf = zeroStart(meanWf)
    %ZEROSTART subtract the first row from all rows
    shape = size(meanWf);
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape(1), []);
    end

    meanWf = bsxfun(@minus, meanWf, meanWf(1, :));
    if numel(shape) ~= 2
        meanWf = reshape(meanWf, shape);
    end
end

function mi = rankorder_mr_(meanWf)
    if isrow(meanWf)
        meanWf = meanWf';
    end

    shape = size(meanWf);
    if numel(shape) == 3
        meanWf = meanWf(:); % global order
    end

    mi = zeros(size(meanWf));
    for iCol = 1:size(meanWf, 2)
        col = meanWf(:, iCol);
        pos = find(col > 0);

        if ~isempty(pos)
            [~, argsort] = sort(col(pos), 'ascend');
            mi(pos(argsort), iCol) = 1:numel(pos);
        end

        neg = find(col < 0);
        if ~isempty(neg)
            [~, argsort] = sort(col(neg), 'descend');
            mi(neg(argsort), iCol) = -(1:numel(neg));
        end
    end

    if numel(shape) == 3
        mi = reshape(mi, shape);
    end
end

function [viSite, viSite2, viSite3] = getSecondaryPeaks(hClust)
    mrMin_clu = squeeze_(min(hClust.tmrWav_spk_clu) - hClust.tmrWav_spk_clu(1,:,:));

    [~, viSite] = min(mrMin_clu);

    mrMin_clu(sub2ind(size(mrMin_clu), viSite, 1:numel(viSite))) = 0;
    [~, viSite2] = min(mrMin_clu);

    mrMin_clu(sub2ind(size(mrMin_clu), viSite2, 1:numel(viSite2))) = 0;
    [~, viSite3] = min(mrMin_clu);
end