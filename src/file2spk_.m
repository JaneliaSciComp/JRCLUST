%--------------------------------------------------------------------------
function S0 = file2spk_(P, viTime_spk0, viSite_spk0)
    % function [tnWav_raw, tnWav_spk, trFet_spk, S0] = file2spk_(P, viTime_spk0, viSite_spk0)
    % file loading routine. keep spike waveform (tnWav_spk) in memory
    % assume that the file is chan x time format
    % usage:
    % [tnWav_raw, tnWav_spk, S0] = file2spk_(P)
    %
    % [tnWav_raw, tnWav_spk, S0] = file2spk_(P, viTime_spk, viSite_spk)
    %   construct spike waveforms from previous time markers
    % 6/29/17 JJJ: Added support for the matched filter

    if nargin < 2
        viTime_spk0 = [];
    end
    if nargin < 3
        viSite_spk0 = [];
    end
    S0 = [];
    % [tnWav_raw, tnWav_spk, trFet_spk, S0] = deal([]);

    viTime_spk0 = viTime_spk0(:);
    viSite_spk0 = viSite_spk0(:);

    if isempty(P.csFile_merge)
        if ~exist_file_(P.vcFile), P.vcFile = subsDir_(P.vcFile, P.vcFile_prm); end
        csFile = {P.vcFile};
    else
        csFile = filter_files_(P.csFile_merge);
        if isempty(csFile)
            P.csFile_merge = subsDir_(P.csFile_merge, P.vcFile_prm);
            csFile = filter_files_(P.csFile_merge);
        end
    end
    if ~isempty(get_(P, 'vcFile_thresh'))
        try
            S_thresh = load(P.vcFile_thresh);
            vnThresh_site = S_thresh.vnThresh_site;
            set0_(vnThresh_site);
            fprintf('Loaded %s\n', P.vcFile_thresh);
        catch
            disperr_('vcFile_thresh load error');
        end
    end
    if isempty(csFile), error('No binary files found.'); end
    % [tnWav_raw, tnWav_spk, trFet_spk, miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site] = deal({});
    [miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site] = deal({});
    viT_offset_file = zeros(size(csFile));
    nFiles = numel(csFile);
    [nSamples1, nLoads] = deal(0); % initialize the counter
    [vrFilt_spk, mrPv_global] = deal([]); % reset the template
    set0_(mrPv_global, vrFilt_spk); % reeset mrPv_global and force it to recompute
    write_spk_(P.vcFile_prm);
    for iFile=1:nFiles
        fprintf('File %d/%d: detecting spikes from %s\n', iFile, nFiles, csFile{iFile});
        t1 = tic;
        [fid1, nBytes_file1] = fopen_(csFile{iFile}, 'r');

        nBytes_file1 = file_trim_(fid1, nBytes_file1, P);

        [nLoad1, nSamples_load1, nSamples_last1] = plan_load_(nBytes_file1, P);
        %         nSamples1 = 0; %accumulated sample offset
        viT_offset_file(iFile) = nSamples1;
        mnWav11_pre = [];
        for iLoad1 = 1:nLoad1
            fprintf('Processing %d/%d of file %d/%d...\n', iLoad1, nLoad1, iFile, nFiles);
            nSamples11 = ifeq_(iLoad1 == nLoad1, nSamples_last1, nSamples_load1);
            fprintf('\tLoading from file...'); t_load_ = tic;
            [mnWav11, vrWav_mean11] = load_file_(fid1, nSamples11, P);
            fprintf('took %0.1fs\n', toc(t_load_));
            if iLoad1 < nLoad1
                mnWav11_post = load_file_preview_(fid1, P);
            else
                mnWav11_post = [];
            end
            [viTime_spk11, viSite_spk11] = filter_spikes_(viTime_spk0, viSite_spk0, nSamples1 + [1, nSamples11]);
            [tnWav_raw_, tnWav_spk_, trFet_spk_, miSite_spk{end+1}, viTime_spk{end+1}, vrAmp_spk{end+1}, vnThresh_site{end+1}, P.fGpu] ...
            = wav2spk_(mnWav11, vrWav_mean11, P, viTime_spk11, viSite_spk11, mnWav11_pre, mnWav11_post);
            write_spk_(tnWav_raw_, tnWav_spk_, trFet_spk_);
            viTime_spk{end} = viTime_spk{end} + nSamples1;
            nSamples1 = nSamples1 + nSamples11;
            if iLoad1 < nLoad1, mnWav11_pre = mnWav11(end-P.nPad_filt+1:end, :); end
            clear mnWav11 vrWav_mean11;
            nLoads = nLoads + 1;
        end %for
        fclose(fid1);
        t_dur1 = toc(t1);
        t_rec1 = (nBytes_file1 / bytesPerSample_(P.vcDataType) / P.nChans) / P.sRateHz;
        fprintf('File %d/%d took %0.1fs (%0.1f MB, %0.1f MB/s, x%0.1f realtime)\n', ...
        iFile, nFiles, ...
        t_dur1, nBytes_file1/1e6, nBytes_file1/t_dur1/1e6, t_rec1/t_dur1);
    end %for
    write_spk_();

    [miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site] = ...
    multifun_(@(x)cat(1, x{:}), miSite_spk, viTime_spk, vrAmp_spk, vnThresh_site);
    vrThresh_site = mean(single(vnThresh_site),1);
    viSite_spk = miSite_spk(:,1);
    if size(miSite_spk,2) >= 2
        viSite2_spk = miSite_spk(:,2);
    else
        viSite2_spk = [];
    end

    % set S0
    [dimm_raw, dimm_spk, dimm_fet] = deal(size(tnWav_raw_), size(tnWav_spk_), size(trFet_spk_));
    [dimm_raw(3), dimm_spk(3), dimm_fet(3)] = deal(numel(viTime_spk));
    nSites = numel(P.viSite2Chan);
    cviSpk_site = arrayfun(@(iSite)find(miSite_spk(:,1) == iSite), 1:nSites, 'UniformOutput', 0);
    if size(miSite_spk,2) >= 2
        cviSpk2_site = arrayfun(@(iSite)find(miSite_spk(:,2) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk2_site = cell(1, nSites);
    end
    if size(miSite_spk,2) >= 3
        cviSpk3_site = arrayfun(@(iSite)find(miSite_spk(:,3) == iSite), 1:nSites, 'UniformOutput', 0);
    else
        cviSpk3_site = [];
    end
    [mrPv_global, vrD_global] = get0_('mrPv_global', 'vrD_global');
    S0 = makeStruct_(P, viSite_spk, viSite2_spk, viTime_spk, vrAmp_spk, vrThresh_site, dimm_spk, ...
    cviSpk_site, cviSpk2_site, cviSpk3_site, dimm_raw, viT_offset_file, dimm_fet, nLoads, ...
    mrPv_global, vrFilt_spk, vrD_global);
end %func
