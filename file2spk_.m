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
    [spikePrSecSites, spikeTimes, vrAmp_spk, siteThresholds] = deal({});
    fileSampleOffsets = zeros(size(filenames));
    [nSamples1, nLoads] = deal(0);
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
        [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes, P);

        fileSampleOffsets(iFile) = nSamples1;
        mnWav11_pre = [];
        for iLoad1 = 1:nLoad1
            fprintf('Processing %d/%d of file %d/%d...\n', iLoad1, nLoad1, iFile, nFiles);
            nSamples11 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);

            fprintf('\tLoading from file...');
            tLoad = tic;
            [mnWav11, vrWav_mean11] = load_file_(fid, nSamples11, P);
            fprintf('took %0.1fs\n', toc(tLoad));

            if iLoad1 < nLoad1
                mnWav11_post = load_file_preview_(fid, P);
            else
                mnWav11_post = [];
            end

            % if given spike times and sites, get the subset of these lying between nSamples1 + 1 and nSamples1 + nSamples11
            [spikeTimes11, spikeSites11] = getSpikesInInterval(spikeTimes0, spikeSites0, nSamples1 + [1, nSamples11]);

            [spikeTraces_, spikeWaveforms_, spikeFeatures_, spikePrSecSites{end+1}, spikeTimes{end+1}, vrAmp_spk{end+1}, siteThresholds{end+1}, P.useGPU] ...
                = wav2spk_(mnWav11, vrWav_mean11, P, spikeTimes11, spikeSites11, mnWav11_pre, mnWav11_post);

            fwrite_(fidTraces, spikeTraces_);
            if strcmp(get_set_(P, 'algorithm', 'JRCLUST'), 'JRCLUST')
                fwrite_(fidWaveforms, spikeWaveforms_);
                fwrite_(fidFeatures, spikeFeatures_);
            end

            spikeTimes{end} = spikeTimes{end} + nSamples1;
            nSamples1 = nSamples1 + nSamples11;

            if iLoad1 < nLoad1
                mnWav11_pre = mnWav11(end-P.nPad_filt+1:end, :);
            end

            clear mnWav11 vrWav_mean11;
            nLoads = nLoads + 1;
        end %for
        fclose(fid);

        t_dur1 = toc(tFile);
        t_rec1 = (nBytes / bytesPerSample_(P.dataType) / P.nChans) / P.sRateHz;

        fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
            iFile, nFiles, t_dur1, nBytes/1e6, nBytes/(t_dur1*1e6), t_rec1/t_dur1);
    end %for

    % close data files
    fclose(fidTraces);
    fclose(fidWaveforms);
    fclose(fidFeatures)

    [spikePrSecSites, spikeTimes, vrAmp_spk, siteThresholds] = ...
        multifun_(@(x) cat(1, x{:}), spikePrSecSites, spikeTimes, vrAmp_spk, siteThresholds);
    vrThresh_site = mean(single(siteThresholds),1);

    spikeSites = spikePrSecSites(:, 1);
    if size(spikePrSecSites, 2) >= 2
        spikeSecondarySites = spikePrSecSites(:,2);
    else
        spikeSecondarySites = [];
    end

    % set S0
    [traceDims, waveformDims, featureDims] = deal(size(spikeTraces_), size(spikeWaveforms_), size(spikeFeatures_));
    [traceDims(3), waveformDims(3), featureDims(3)] = deal(numel(spikeTimes));
    nSites = numel(P.chanMap);
    cviSpk_site = arrayfun(@(iSite)find(spikePrSecSites(:,1) == iSite), 1:nSites, 'UniformOutput', 0);
    if size(spikePrSecSites, 2) >= 2
        cviSpk2_site = arrayfun(@(iSite)find(spikePrSecSites(:,2) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk2_site = cell(1, nSites);
    end
    if size(spikePrSecSites,2) >= 3
        cviSpk3_site = arrayfun(@(iSite)find(spikePrSecSites(:,3) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk3_site = [];
    end
    [mrPv_global, vrD_global] = get0_('mrPv_global', 'vrD_global');

    % save everything
    S0 = makeStruct_(P, spikeSites, spikeSecondarySites, spikeTimes, vrAmp_spk, vrThresh_site, waveformDims, ...
        cviSpk_site, cviSpk2_site, cviSpk3_site, traceDims, fileSampleOffsets, featureDims, nLoads, ...
        mrPv_global, vrFilt_spk, vrD_global);
end %func
