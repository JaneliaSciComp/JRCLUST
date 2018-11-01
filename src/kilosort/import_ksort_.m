%--------------------------------------------------------------------------
function import_ksort_(vcFile_prm)
    %IMPORT_KSORT_ import a Kilosort(2) session from rez.mat

    % short-circuit standard param loading
    P = file2struct_(vcFile_prm);
    vcName_session = strrep(vcFile_prm, '.prm', '');

    if isempty(P)
        P = struct();
    end

    if ~isfield(P, 'template_file')
        P.template_file = '';
    elseif ~isempty(P.template_file)
        assert_(exist_file_(P.template_file), sprintf('template file does not exist: %s', P.template_file));
        P = struct_merge_(file2struct_(P.template_file), P);
    end
    P.vcFile_prm = vcFile_prm;

    if ~isfield(P, 'vcFile_rez') || isempty(P.vcFile_rez)
        error('define vcFile_rez in your .prm file');
    end

    try
        rez = load(P.vcFile_rez, '-mat', 'rez');
    catch err
        if startsWith(err.message, 'Unable to read MAT-file')
            error('rez file must be a MAT file');
        else
            rethrow(err);
        end
    end

    if isempty(rez)
        error('rez is empty');
    end

    if isfield(rez, 'rez')
        rez = rez.rez;
    end

    if isfield(rez, 'Wphy')
        rez.W = rez.Wphy;
    end

    % these matrices are often (always?) stored as GPU arrays
    rez.W = gather(rez.W);
    rez.U = gather(rez.U);
    rez.mu = gather(rez.mu);

    % import spike location, time and cluster assignment
    viTime_spk = rez.st3(:, 1); % spike times
    viTemplate_spk = rez.st3(:, 2); % spike templates

    if size(rez.st3, 2) == 5
        viClu = 1 + rez.st3(:, 5); % post-merging result
    else
        viClu = viTemplate_spk; % template/cluster
    end

    [clusters, ~, indices] = unique(viClu);
    % separate out the good ones from the junk ones
    viClu_pos = clusters(clusters > 0);
    viClu_junk = setdiff(clusters, viClu_pos);
    viClu_new = [viClu_junk' 1:numel(viClu_pos)]';
    viClu = viClu_new(indices);

    nTemplates = size(rez.simScore, 1);
    nClu = numel(viClu_pos);

    cviTemplate_clu = arrayfun(@(iClu) unique(viTemplate_spk(viClu == iClu)), 1:nClu, 'UniformOutput', 0);

    % compute templates
    nt0 = size(rez.W, 1);
    U = rez.U;
    W = rez.W;

    Nfilt = size(W, 2);
    Nchan = rez.ops.Nchan;

    trTemplates = zeros(Nchan, nt0, Nfilt, 'single');
    for iNN = 1:size(trTemplates,3)
       trTemplates(:,:,iNN) = squeeze(U(:,iNN,:)) * squeeze(W(:,iNN,:))';
    end
    trTemplates = -abs(permute(trTemplates, [3 2 1])); % nTemplates x nSamples x nChannels

    % compute the weighted average template for a given cluster and pick its min site
    viSite_clu = zeros('like', clusters);
    for iClu = 1:nClu
        viTemp_clu = cviTemplate_clu{iClu}; % templates for this cluster

        mrTemp_mean = squeeze(trTemplates(viTemp_clu, :, :));
        if numel(viTemp_clu) > 1
            freqs = histcounts(viTemplate_spk(viClu == iClu), numel(viTemp_clu));
            weights = freqs/sum(freqs); % weight sum of templates by frequency of occurrence in this cluster
            t = zeros(size(mrTemp_mean, 2), size(mrTemp_mean, 3));
            for iWeight = numel(weights)
                t = t + weights(iWeight)*squeeze(mrTemp_mean(iWeight, :, :));
            end

            mrTemp_mean = t;
        end

        sampleMin = min(mrTemp_mean, [], 1);
        [~, viSite_clu(iClu)] = min(sampleMin); % cluster location
    end
    viSite_spk = viSite_clu(viClu);

    % construct P from scratch
    P.fImportKsort = 1;
    P.vcFet_show = 'kilosort';
    P.vcFile = get_set_(P, 'vcFile', rez.ops.fbinary);
    P.csFile_merge = ''; % kilosort takes concatenated files
    P.sRateHz = rez.ops.fs;

    % probe
    P.viSite2Chan = rez.ops.chanMap;
    P.mrSiteXY = [rez.xcoords(:) rez.ycoords(:)];
    P.nChans = rez.ops.NchanTOT;
    P.viShank_site = rez.ops.kcoords(:)';
    P.vrSiteHW = get_set_(P, 'vrSiteHW', [12 12]); % TODO: address this
    if isfield(rez, 'connected')
        P.viSite2Chan = P.viSite2Chan(rez.connected);
        P.mrSiteXY = P.mrSiteXY(rez.connected(:), :);
    end
    P.viChan_aux = setdiff(1:P.nChans, 1:max(P.viSite2Chan));

    % loading settings
    P.MAX_LOAD_SEC = get_set_(P, 'MAX_LOAD_SEC', []);

    % filter settings
    P.vcFilter = get_set_(P, 'vcFilter', 'bandpass');
    P.freqLim = get_set_(P, 'freqLim', [500, .475 * P.sRateHz]);

    P.fft_thresh = get_set_(P, 'fft_thresh', 0);
    P.vcCommonRef = get_set_(P, 'vcCommonRef', 'mean');
    P.blank_thresh = get_set_(P, 'blank_thresh', 0);
    P.vcDataType = get_set_(P, 'vcDataType', 'int16'); % KS default

    P.tlim_load = get_set_(P, 'tlim_load', []);
    P.fTranspose_bin = get_set_(P, 'fTranspose_bin', 1);
    P.maxSite = get_set_(P, 'maxSite', 6.5); % TODO: allow user to set
    P.nSites_ref = get_set_(P, 'nSites_ref', 0); % TODO: address

    P.spkLim_ms = get_set_(P, 'spkLim_ms', [-.25 .75]); % JRC default
    P.spkLim = get_set_(P, 'spkLim', round(P.spkLim_ms/1000*P.sRateHz)); % JRC default
    P.spkLim_raw = get_set_(P, 'spkLim_raw', ceil(size(trTemplates, 2)/2) * [-1, 1]); % TODO: address

    P.corrLim = get_set_(P, 'corrLim', [.75 1]);
    P.fDrift_merge = get_set_(P, 'fDrift_merge', 0); % do not attempt drift correction
    P.nPaddingSamples = get_set_(P, 'nPaddingSamples', 100); % JRC default
    P.qqFactor = get_set_(P, 'qqFactor', 5); % JRC default
    P.nPcPerChan = get_set_(P, 'nPcPerChan', 1); % JRC default
    P.nTime_clu = get_set_(P, 'nTime_clu', 1); % spikes detected and clustered over entire time series
    P.uV_per_bit = get_set_(P, 'uV_per_bit', 0.305176); % JRC default
    P.spkRefrac_ms = get_set_(P, 'spkRefrac_ms', .25); % JRC default

    P.viSiteZero = get_set_(P, 'viSiteZero', []);
    P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site);
    P.fGPU = double(gpuDeviceCount() > 0);

    P.vcFet = get_set_(P, 'vcFet', 'gpca'); % JRC default

    S0 = file2spk_(P, int32(viTime_spk), int32(viSite_spk));
    P = save_prb_([vcName_session '-probe.mat'], P);
    S0.P = P;
    S0.rez = rez;

    set(0, 'UserData', S0);

    global trFet_spk
    trFet_spk = get_spkfet_(P);

    % construct S_clu from scratch
    S_clu = struct();
    S_clu.nTemplates = nTemplates;
    S_clu.viTemplate_spk = viTemplate_spk;
    S_clu.cviTemplate_clu = cviTemplate_clu;
    S_clu.trTemplates = trTemplates;

    S_clu.nClu = nClu;
    S_clu.viClu = viClu;
    S_clu.viClu_auto = viTemplate_spk;

    S_clu.csNote_clu = cell(nClu, 1);
    S_clu.viSite_clu = viSite_clu;

    S_clu.cviSpk_clu = arrayfun(@(iClu) find(viClu == iClu), 1:nClu, 'UniformOutput', 0);

    S_clu = S_clu_refresh_(S_clu, 0); % don't remove empty
    S_clu = S_clu_sort_(S_clu, 'viSite_clu');
    S_clu = S_clu_update_wav_(S_clu);
    S_clu.P = P;
    S_clu = S_clu_position_(S_clu);
    S_clu = S_clu_update_(S_clu, 1:nClu, P);
    S_clu = sim_score_(S_clu);

    assert_(S_clu_valid_(S_clu), 'Import failed: inconsistent clusters.');
    S0.S_clu = S_clu;

    % Save
    set(0, 'UserData', S0);

    export_prm_(P.vcFile_prm, [], 0);
    save0_([vcName_session '_jrc.mat']);
end % function
