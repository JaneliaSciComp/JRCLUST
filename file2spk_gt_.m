%--------------------------------------------------------------------------
function [spikeWaveforms, vrVrms_site] = file2spk_gt_(P, spikeTimes0)
    % file loading routine. keep spike waveform (spikeWaveforms) in memory
    % assume that the file is chan x time format
    % usage:
    % [spikeTraces, spikeWaveforms, S0] = file2spk_(P)
    %
    % [spikeTraces, spikeWaveforms, S0] = file2spk_(P, spikeTimes, spikeSites)
    %   construct spike waveforms from previous time markers
    % 6/29/17 JJJ: Added support for the matched filter
    P.fft_thresh = 0; %disable for GT
    [spikeWaveforms, siteThresholds] = deal({});
    nSamples1 = 0;
    [fid1, nBytes_file1] = fopenInfo(P.vcFile, 'r');
    nBytes_file1 = file_trim_(fid1, nBytes_file1, P);
    [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file1, P);
    t_dur1 = tic;
    mnWav11_pre = [];
    for iLoad1 = 1:nLoad1
        fprintf('\tProcessing %d/%d...\n', iLoad1, nLoad1);
        nSamples11 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);
        [mnWav11, vrWav_mean11] = load_file_(fid1, nSamples11, P);
        if iLoad1 < nLoad1
            mnWav11_post = load_file_preview_(fid1, P);
        else
            mnWav11_post = [];
        end
        [spikeTimes11] = filter_spikes_(spikeTimes0, [], nSamples1 + [1, nSamples11]);
        [spikeWaveforms{end+1}, siteThresholds{end+1}] = wav2spk_gt_(mnWav11, P, spikeTimes11, mnWav11_pre, mnWav11_post);
        if iLoad1 < nLoad1, mnWav11_pre = mnWav11(end-P.nPad_filt+1:end, :); end
        nSamples1 = nSamples1 + nSamples11;
        clear mnWav11 vrWav_mean11;
    end %for
    fclose(fid1);
    t_dur1 = toc(t_dur1);
    t_rec1 = (nBytes_file1 / bytesPerSample_(P.dataType) / P.nChans) / P.sRateHz;
    fprintf('took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
    t_dur1, nBytes_file1/1e6, nBytes_file1/t_dur1/1e6, t_rec1/t_dur1);

    % spikeTraces = cat(3, spikeTraces{:});
    spikeWaveforms = cat(3, spikeWaveforms{:});
    [siteThresholds] = multifun_(@(x)cat(1, x{:}), siteThresholds);
    vrVrms_site = mean(single(siteThresholds),1) / P.qqFactor;

end %func
