%--------------------------------------------------------------------------
% 2017/12/5 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_(S_clu, P, viClu_update)
    % version selector
    fUse_old = 1;

    spkLim_raw_factor = get_set_(P, 'spkLim_raw_factor', 2);
    spkLim_factor_merge = get_set_(P, 'spkLim_factor_merge', 1);
    if nargin<3, viClu_update = []; end

    if spkLim_factor_merge < spkLim_raw_factor && ~fUse_old
        mrWavCor = S_clu_wavcor_2_(S_clu, P, viClu_update); % newer format v3.1.8
    else
        mrWavCor = S_clu_wavcor_1_(S_clu, P, viClu_update); % older format
    end
end %func

%% local functions
%--------------------------------------------------------------------------
% 6/29/17 JJJ: distance-based neighboring unit selection
function mrWavCor = S_clu_wavcor_1_(S_clu, P, viClu_update)
    % symmetric matrix and common basis comparison only

    fUsePeak2 = 0;
    nInterp_merge = get_set_(P, 'nInterp_merge', 1); % set to 1 to disable
    fDrift_merge = get_set_(P, 'fDrift_merge', 0);
    P.fGpu = 0;
    nShift = ceil(P.spkRefrac_ms / 1000 * P.sRateHz * nInterp_merge); % +/-n number of samples to compare time shift
    % nShift = 0;
    fWaveform_raw = get_set_(P, 'fWavRaw_merge', 1); % revert TW: get_set_(P, 'fWavRaw_merge', 0);

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
    % tmrWav_clu = gpuArray_(tmrWav_clu, P.fGpu);
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
    mrWavCor = gpuArray_(zeros(nClu), P.fGpu);
    nSites_spk = P.maxSite*2+1-P.nSites_ref;
    nT = size(tmrWav_clu, 1);
    [cviShift1, cviShift2] = shift_range_(gpuArray_(int32(nT), P.fGpu), nShift);
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
    nClu = numel(viSite_clu);
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

    vrWavCor2 = zeros(nClu, 1, 'single');
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

    nClu = S_clu.nClu;
    mrWavCor = gpuArray_(zeros(nClu), P.fGpu);
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

%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, viSite_clu, iClu2, cell_5args)
    [P, vlClu_update, mrWavCor0, cviShift1, cviShift2] = deal(cell_5args{:});
    nClu = numel(viSite_clu);
    iSite_clu2 = viSite_clu(iClu2);
    if iSite_clu2==0 || isnan(iSite_clu2), vrWavCor2 = []; return; end
    viSite2 = P.miSites(:,iSite_clu2);
    maxDist_site_um = get_set_(P, 'maxDist_site_merge_um', 35);
    viClu1 = find(ismember(viSite_clu, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));

    vrWavCor2 = zeros(nClu, 1, 'single');
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
            mr1_ = reshape(meanSubt_(tr1(cviShift1{iShift},:,:)), [], nDrift);
            mr2_ = reshape(meanSubt_(tr2(cviShift2{iShift},:,:)), [], nDrift);

            mr12_ = (mr1_' * mr2_) ./ sqrt(sum(mr1_.^2)' * sum(mr2_.^2));
            vrCor(iShift) = max(mr12_(:));
        end
        maxCor = max(vrCor);
    end
end
