%--------------------------------------------------------------------------
function S0 = file2spk_(P, spikeTimes0, spikeSites0)
    % function [spikeTraces, spikeWaveforms, spikeFeatures, S0] = file2spk_(P, spikeTimes0, spikeSites0)
    % file loading routine. keep spike waveform (spikeWaveforms) in memory
    % assume that the file is chan x time format
    % usage:
    % [spikeTraces, spikeWaveforms, S0] = file2spk_(P)
    %
    % [spikeTraces, spikeWaveforms, S0] = file2spk_(P, spikeTimes, spikeSites)
    %   construct spike waveforms from previous time markers
    % 6/29/17 JJJ: Added support for the matched filter

    if nargin < 2
        spikeTimes0 = [];
    end

    if nargin < 3
        spikeSites0 = [];
    end

    S0 = [];

    spikeTimes0 = spikeTimes0(:);
    spikeSites0 = spikeSites0(:);

    % regularize list of files to load
    if ~isfield(P, 'multiFilenames') || isempty(P.multiFilenames)
        if ~fileExists(P.vcFile)
            P.vcFile = replaceDir(P.vcFile, P.paramFile);
        end

        filenames = {P.vcFile};
    else
        filenames = filter_files_(P.multiFilenames);
        if isempty(filenames)
            P.multiFilenames = replaceDir(P.multiFilenames, P.paramFile);
            filenames = filter_files_(P.multiFilenames);
        end
    end

    % load site-wise thresholds, if applicable
    if ~isempty(get_(P, 'thresholdFile'))
        try
            S_thresh = load(P.thresholdFile);
            siteThresholds = S_thresh.siteThresholds;
            setUserData(siteThresholds);
            fprintf('Loaded %s\n', P.thresholdFile);
        catch
            disperr_('thresholdFile load error');
        end
    end

    if isempty(filenames)
        error('No binary files found.');
    end

    nFiles = numel(filenames);

    % initialize lists and counters
    [spikePrimarySecondarySites, spikeTimes, vrAmp_spk, siteThresholds] = deal({});
    fileSampleOffsets = zeros(size(filenames));
    [offset, nLoads] = deal(0);
    [vrFilt_spk, mrPv_global] = deal([]);

    setUserData(mrPv_global, vrFilt_spk); % reset mrPv_global and force it to recompute

    fidTraces = fopen(strrep(P.paramFile, '.prm', '_traces.bin'), 'W');
    fidWaveforms = fopen(strrep(P.paramFile, '.prm', '_waveforms.bin'), 'W');
    fidFeatures = fopen(strrep(P.paramFile, '.prm', '_features.bin'), 'W');

    for iFile = 1:nFiles
        fprintf('File %d/%d: detecting spikes from %s\n', iFile, nFiles, filenames{iFile});
        tFile = tic;

        [fid, nBytes] = fopenInfo(filenames{iFile}, 'r');
        nBytes = file_trim_(fid, nBytes, P);
        [nFileLoads, nSamples_load1, nSamples_last1] = plan_load_(nBytes, P);

        fileSampleOffsets(iFile) = offset;
        prePadding = [];

        for iLoad = 1:nFileLoads
            fprintf('Processing %d/%d of file %d/%d...\n', iLoad, nFileLoads, iFile, nFiles);

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

            if iLoad < nFileLoads && P.nPaddingSamples > 0
                postPadding = loadFilePreview(fid, P);
            else
                postPadding = [];
            end

            [spikeTraces_, spikeWaveforms_, spikeFeatures_, spikePrimarySecondarySites{end+1}, spikeTimes{end+1}, vrAmp_spk{end+1}, siteThresholds{end+1}, P.useGPU] ...
                = wav2spk_(loadSamples, channelMeans, P, loadSpikeTimes, loadSpikeSites, prePadding, postPadding);

            fwrite_(fidTraces, spikeTraces_);
            if get_set_(P, 'fImportKilosort', 0)
                fwrite_(fidWaveforms, spikeWaveforms_);
                fwrite_(fidFeatures, spikeFeatures_);
            end

            spikeTimes{end} = spikeTimes{end} + offset;
            offset = offset + nLoadSamples;

            if iLoad < nFileLoads
                prePadding = loadSamples(end - P.nPaddingSamples + 1:end, :);
            end

            clear loadSamples channelMeans;
            nLoads = nLoads + 1;
        end %for
        fclose(fid);

        t_dur1 = toc(tFile);
        t_rec1 = (nBytes / bytesPerSample_(P.dataType) / P.nChans) / P.sampleRateHz;

        fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
            iFile, nFiles, t_dur1, nBytes*1e-6, nBytes/(t_dur1*1e6), t_rec1/t_dur1);
    end %for

    % close data files
    fclose(fidTraces);
    fclose(fidWaveforms);
    fclose(fidFeatures)

    [spikePrimarySecondarySites, spikeTimes, vrAmp_spk, siteThresholds] = ...
        multifun_(@(x) cat(1, x{:}), spikePrimarySecondarySites, spikeTimes, vrAmp_spk, siteThresholds);
    vrThresh_site = mean(single(siteThresholds),1);

    spikeSites = spikePrimarySecondarySites(:, 1);
    if size(spikePrimarySecondarySites, 2) >= 2
        spikeSecondarySites = spikePrimarySecondarySites(:,2);
    else
        spikeSecondarySites = [];
    end

    % set S0
    [traceDims, waveformDims, featureDims] = deal(size(spikeTraces_), size(spikeWaveforms_), size(spikeFeatures_));
    [traceDims(3), waveformDims(3), featureDims(3)] = deal(numel(spikeTimes));
    nSites = numel(P.chanMap);
    cviSpk_site = arrayfun(@(iSite)find(spikePrimarySecondarySites(:,1) == iSite), 1:nSites, 'UniformOutput', 0);
    if size(spikePrimarySecondarySites, 2) >= 2
        cviSpk2_site = arrayfun(@(iSite)find(spikePrimarySecondarySites(:,2) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk2_site = cell(1, nSites);
    end
    if size(spikePrimarySecondarySites,2) >= 3
        cviSpk3_site = arrayfun(@(iSite)find(spikePrimarySecondarySites(:,3) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk3_site = [];
    end
    [mrPv_global, vrD_global] = get0_('mrPv_global', 'vrD_global');

    % save everything
    S0 = makeStruct_(P, spikeSites, spikeSecondarySites, spikeTimes, vrAmp_spk, vrThresh_site, waveformDims, ...
        cviSpk_site, cviSpk2_site, cviSpk3_site, traceDims, fileSampleOffsets, featureDims, nLoads, ...
        mrPv_global, vrFilt_spk, vrD_global);
end %func
