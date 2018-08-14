%--------------------------------------------------------------------------
function S0 = kilosort2jrc_(P, spikeTimes0, spikeSites0)
    if nargin < 2
        spikeTimes0 = [];
    end

    if nargin < 3
        spikeSites0 = [];
    end

    S0 = [];

    spikeTimes0 = spikeTimes0(:);
    spikeSites0 = spikeSites0(:);

    filename = P.vcFile;
    if isempty(filename)
        error('No binary file given.');
    end

    % initialize lists and counters
    spikeTimes = {};

    [offset, nLoads] = deal(0);
    [vrFilt_spk, mrPv_global] = deal([]);

    setUserData(mrPv_global, vrFilt_spk); % reset mrPv_global and force it to recompute

    fidTraces = fopen(strrep(P.paramFile, '.prm', '_traces.bin'), 'W');
    fidWaveforms = fopen(strrep(P.paramFile, '.prm', '_waveforms.bin'), 'W');

    fprintf('Extracting spikes from %s\n', filename);
    tFile = tic;

    [fid, nBytes] = fopenInfo(filename, 'r');
    nBytes = file_trim_(fid, nBytes, P);
    [nFileLoads, nSamples_load1, nSamples_last1] = planLoad(nBytes, P);

    for iLoad = 1:nFileLoads
        fprintf('Processing %d/%d of file...\n', iLoad, nFileLoads);

        if iLoad < nFileLoads
            nLoadSamples = nSamples_load1;
        else
            nLoadSamples = nSamples_last1;
        end

        % if given spike times and sites, get the subset of these lying in the interval offset + [1, nLoadSamples]
        [loadSpikeTimes, loadSpikeSites] = getIntervalTimesSites(spikeTimes0, spikeSites0, offset + [1, nLoadSamples]);

        fprintf('\tLoading from file...');
        tLoad = tic;
        [loadSamples, channelMeans] = load_file_(fid, nLoadSamples, P);
        fprintf('took %0.1fs\n', toc(tLoad));

        % extract traces
        spikeTraces_ = tracesToTensor(loadSamples, loadSpikeSites, loadSpikeTimes, P.spkLim_raw, P);
        fwrite_(fidTraces, spikeTraces_);

        % filter traces (lifted from wav2spk_)
        fprintf('\tFiltering spikes...');
        tFilter = tic;

        loadSamplesCPU = loadSamples; % keep a copy in CPU
        try
            [loadSamples, P.useGPU] = gpuArray_(loadSamples, P.useGPU);
        catch % GPu failure
            P.useGPU = 0;
            loadSamples = loadSamplesCPU;
        end

        clear loadSamplesCPU;

        [loadFilteredSamples, ~] = filt_car_(loadSamples, P);
        fprintf('\ttook %0.1fs\n', toc(tFilter));

        spikeWaveforms_ = tracesToTensor(loadFilteredSamples, loadSpikeSites, loadSpikeTimes, P.spkLim, P);
        fwrite_(fidWaveforms, spikeWaveforms_);

        spikeTimes{end+1} = loadSpikeTimes + offset;
        offset = offset + nLoadSamples;

        clear loadSamples channelMeans;
        nLoads = nLoads + 1;
    end %for

    fclose(fid);

    tFile = toc(tFile);
    tRecording = (nBytes / bytesPerSample_(P.dataType) / P.nChans) / P.sampleRateHz;

    fprintf('File took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
        tFile, nBytes*1e-6, nBytes*1e-6/tFile, tRecording/tFile);

    % close data files
    fclose(fidTraces);
    fclose(fidWaveforms);
    
    spikeTimes = cat(1, spikeTimes{:});

    spikeSites = spikeSites0;
    spikeSecondarySites = [];

    % set S0
    [traceDims, waveformDims] = deal(size(spikeTraces_), size(spikeWaveforms_));
    [traceDims(3), waveformDims(3)] = deal(numel(spikeTimes));
    nSites = numel(P.chanMap);

    cviSpk_site = arrayfun(@(iSite) find(spikeSites(:, 1) == iSite), 1:nSites, 'UniformOutput', 0);
    cviSpk2_site = cell(1, nSites);
    cviSpk3_site = [];

    [mrPv_global, vrD_global] = get0_('mrPv_global', 'vrD_global');

    % save everything
    S0 = makeStruct_(P, spikeSites, spikeSecondarySites, spikeTimes, waveformDims, ... % vrAmp_spk, vrThresh_site
        cviSpk_site, cviSpk2_site, cviSpk3_site, traceDims, nLoads, mrPv_global, vrFilt_spk, vrD_global);
end %func
