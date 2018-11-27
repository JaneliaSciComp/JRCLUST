%--------------------------------------------------------------------------
% 10/26/17 JJJ: Created
function import_h5_(vcFile_h5)
    % Brian Allen h5 data
    if nargin<1, vcFile_h5 = ''; end
    if isempty(vcFile_h5), vcFile_h5 = 'E:\BrianAllen\915_18\915_18_1\915_18_1.h5'; end
    [vcDir_, vcFile_, vcExt_] = fileparts(vcFile_h5);
    if isempty(vcExt_)
        vcFile_h5 = fullfile(vcFile_h5, [vcFile_, '.h5']);
    end

    P = struct('vcFile', strrep(vcFile_h5, '.h5', '.bin'), 'qqFactor', 5, ...
    'maxDist_site_um', 50, 'maxDist_site_spk_um', 70, 'uV_per_bit', .195, ...
    'max_bursting_index', 3, 'nTime_clu', 4, 'fft_thresh', 0); % set to [] to disable
    S_gt = struct();
    S_gt.probe_layout           = h5readatt(vcFile_h5, '/','probelayout');
    try P.viChanZero           = h5readatt(vcFile_h5, '/','badchannels'); catch, end
    S_gt.sRateHz_gt = h5readatt(vcFile_h5, '/', 'abfsamplerate');
    S_gt.patchtype = h5readatt(vcFile_h5, '/', 'patchtype');
    try  S_gt.padmaptextname = h5readatt(vcFile_h5, '/', 'padmaptextname'); catch, end
    try S_gt.padpitch = h5readatt(vcFile_h5, '/', 'padpitch'); catch, end

    P.sRateHz        = h5readatt(vcFile_h5, '/', 'MEAsamplerate');
    P.nChans                = S_gt.probe_layout(1) * S_gt.probe_layout(2);
    % P = h5readatt_(vcFile_h5, {'patchtype', 'padmaptextname', 'patchtype', 'badchannels'});

    [vcDir, vcFile, vcExt] = fileparts(vcFile_h5);
    vcFile_raw = fullfile(vcDir, 'Recordings', [vcFile, '_raw.h5']);
    vcFile_filtered = fullfile(vcDir, 'Recordings', [vcFile, '_filtered.h5']);
    vcFile_spikes = fullfile(vcDir, 'Analyses', [vcFile, '_spikes.h5']);

    % photodiode = h5read(vcFile_raw, '/photodiode');
    % syncMEA = h5read(vcFile_raw, '/syncMEA');
    % filteredMEA = h5read(vcFile_filtered, '/filteredMEA');
    if ~exist_file_(P.vcFile)
        if exist_file_(vcFile_filtered)
            mrWav = h5read(vcFile_filtered, '/filteredMEA');
        elseif exist_file_(vcFile_raw)
            mrWav = h5read(vcFile_raw, '/rawMEA');
        else
            error('no traces found');
        end
        % P.uV_per_bit = min_step_(rawMEA(:,1));
        write_bin_(P.vcFile, int16(jrclust.utils.meanSubtract(mrWav) / P.uV_per_bit)');
    end


    % Create GT
    if exist_file_(vcFile_spikes)
        S_gt.viTime_all = ceil(h5read(vcFile_spikes, '/derivspiketimes') * S_gt.sRateHz_gt);
        S_gt.vrBI_all = h5read(vcFile_spikes, '/burstindex');
        if ~isempty(get_(P, 'max_bursting_index'))
            viTime_gt = S_gt.viTime_all(S_gt.vrBI_all < P.max_bursting_index); % non-bursting only
        else
            viTime_gt = S_gt.viTime_all;
        end
        %     rawPipette = h5read(vcFile_filtered, '/rawPipette');
    elseif exist_file_(vcFile_raw)
        rawPipette = h5read(vcFile_raw, '/rawPipette');
        rawPipette1 = ndiff_(rawPipette, 2);
        [viTime_gt, vrAmp_gt, thresh_gt] = spikeDetectSingle_fast_(-rawPipette1, struct('qqFactor', 10));
        %     figure; hold on; plot(rawPipette1); plot(S_gt.viTime_all, rawPipette1(S_gt.viTime_all), 'ro');
    else
        error('no spike time info');
    end

    S_gt.viTime = viTime_gt;
    S_gt.viClu = ones(size(viTime_gt));
    jrclust.utils.saveStruct(S_gt, strrep(P.vcFile, '.bin', '_gt.mat'));


    % Create prm file
    P.probe_file = sprintf('boyden%d.prb', P.nChans);
    P.vcFile_prm = [strrep(P.vcFile, '.bin', '_'), strrep(P.probe_file, '.prb', '.prm')];
    copyfile(default_prm_path_(), P.vcFile_prm, 'f');
    edit_prm_file_(P, P.vcFile_prm);
    assignWorkspace_(P, S_gt);
    fprintf('Created .prm file: %s\n', P.vcFile_prm);
    edit(P.vcFile_prm);
    jrc('setprm', P.vcFile_prm); % set the currently working prm file
end %func
