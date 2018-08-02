%--------------------------------------------------------------------------
% 10/10/17 JJJ: moved tnWav_spk and tnWav_raw internally
function S_clu = S_clu_wav_(S_clu, viClu_update, fSkipRaw)
    % average cluster waveforms and determine the center
    % only use the centered spikes
    % viClu_copy: waveforms not changing
    if nargin<2, viClu_update = []; end
    if nargin<3, fSkipRaw = 0; end
    % P = S_clu.P;
    % P = get0_('P');
    S0 = get(0, 'UserData');
    P = S0.P;
    % [dimm_spk, dimm_raw] = get0_('dimm_spk', 'dimm_raw');
    P.fMeanSubt = 0;
    fVerbose = isempty(viClu_update);
    if fVerbose, fprintf('Calculating cluster mean waveform.\n\t'); t1 = tic; end
    if isfield(S_clu, 'nClu')
        nClu = S_clu.nClusters;
    else
        nClu = max(S_clu.spikeClusters);
    end
    nSamples = S0.dimm_spk(1);
    nSites = numel(P.chanMap);
    nSites_spk = S0.dimm_spk(2); % n sites per event group (maxSite*2+1);

    % Prepare cluster loop
    trWav_spk_clu = zeros([nSamples, nSites_spk, nClu], 'single');
    % vrFracCenter_clu = zeros(nClu, 1);
    tmrWav_spk_clu = zeros(nSamples, nSites, nClu, 'single');
    if ~fSkipRaw
        nSamples_raw = S0.dimm_raw(1);
        [trWav_raw_clu] = deal(zeros(nSamples_raw, nSites_spk, nClu, 'single'));
        [tmrWav_raw_clu, tmrWav_raw_lo_clu, tmrWav_raw_hi_clu] = deal(zeros(nSamples_raw, nSites, nClu, 'single'));
    else
        [trWav_raw_clu, tmrWav_raw_clu, tmrWav_raw_lo_clu, tmrWav_raw_hi_clu] = deal([]);
    end
    if ~isempty(viClu_update)
        vlClu_update = false(nClu, 1);
        vlClu_update(viClu_update) = 1;
        nClu_pre = size(S_clu.trWav_spk_clu, 3);
        vlClu_update((1:nClu) > nClu_pre) = 1;
        [tmrWav_spk_clu, tmrWav_raw_clu] = deal(S_clu.tmrWav_spk_clu, S_clu.tmrWav_raw_clu);
        [trWav_spk_clu, trWav_raw_clu] = deal(S_clu.trWav_spk_clu, S_clu.trWav_raw_clu);
        [tmrWav_raw_lo_clu, tmrWav_raw_hi_clu] = deal(S_clu.tmrWav_raw_lo_clu, S_clu.tmrWav_raw_hi_clu);
    else
        vlClu_update = true(nClu, 1);
    end

    % Compute spkwav
    tnWav_ = get_spkwav_(P, 0);
    for iClu=1:nClu
        if vlClu_update(iClu)
            [mrWav_clu1, clusterSites1] = clu_wav_(S_clu, tnWav_, iClu, S0);
            if isempty(mrWav_clu1), continue; end
            [tmrWav_spk_clu(:,clusterSites1,iClu), trWav_spk_clu(:,:,iClu)] = ...
            deal(bit2uV_(mrWav_clu1, P));
        end
        if fVerbose, fprintf('.'); end
    end %clu

    % Compute spkraw
    if ~fSkipRaw
        tnWav_ = []; % clear memory
        tnWav_ = get_spkwav_(P, 1);
        for iClu=1:nClu
            if vlClu_update(iClu)
                [mrWav_clu1, clusterSites1, mrWav_lo_clu1, mrWav_hi_clu1] = clu_wav_(S_clu, tnWav_, iClu, S0);
                if isempty(mrWav_clu1), continue; end
                [tmrWav_raw_clu(:,clusterSites1,iClu), trWav_raw_clu(:,:,iClu)] = deal(meanSubt_(mrWav_clu1) * P.uV_per_bit);
                if isempty(mrWav_lo_clu1) || isempty(mrWav_hi_clu1), continue; end
                tmrWav_raw_lo_clu(:,clusterSites1,iClu) = meanSubt_(mrWav_lo_clu1) * P.uV_per_bit;
                tmrWav_raw_hi_clu(:,clusterSites1,iClu) = meanSubt_(mrWav_hi_clu1) * P.uV_per_bit;
            end
            if fVerbose, fprintf('.'); end
        end %clu
    end

    tmrWav_clu = tmrWav_spk_clu; %meanSubt_ after or before?

    % measure waveforms
    [vrVmin_clu, viSite_min_clu] = min(permute(min(trWav_spk_clu),[2,3,1]),[],1);
    vrVmin_clu = abs(vrVmin_clu(:));
    viSite_min_clu = viSite_min_clu(:);

    S_clu = struct_add_(S_clu, vrVmin_clu, viSite_min_clu, ...
    trWav_spk_clu, tmrWav_spk_clu, trWav_raw_clu, tmrWav_raw_clu, tmrWav_clu, ...
    tmrWav_raw_lo_clu, tmrWav_raw_hi_clu);
    if fVerbose, fprintf('\n\ttook %0.1fs\n', toc(t1)); end
end %func
