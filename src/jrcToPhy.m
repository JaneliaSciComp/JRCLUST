function [spikeTimes, spikeClusters] = jrcToPhy(S0, savePath)
    %JRCTOPHY Summary of this function goes here
    %   Detailed explanation goes here
    if nargin < 1
        S0 = get(0, 'UserData');
    end
    if nargin < 2
        savePath = '.';
    end
    
    S0.P.vcFet_show = 'vpp';
    P = S0.P;
    
    S_clu = S0.S_clu;

    tnWav_spk = double(get_spkwav_(P, 0));
    P.viSite2Chan = P.viSite2Chan(:);

    outputs = {'channel_map.npy', 'channel_positions.npy', ...
               'pc_features.npy', 'spike_channels.npy', ...
               'spike_clusters.npy', 'spike_times.npy', ...
               'vpp_features.npy'};

    fs = dir(fullfile(savePath, '*.npy'));
    for i = 1:length(fs)
        fname = fs(i).name;
        % don't delete .npy files which have nothing to do with us
        if find(strcmp(fname, outputs))
            delete(fullfile(savePath, fname));
        end
    end

    if exist(fullfile(savePath, '.phy'), 'dir')
        rmdir(fullfile(savePath, '.phy'), 's');
    end
    
    nSpikes = numel(S0.viTime_spk);
    nSites = numel(P.viSite2Chan);
    
    % data for .npy files
    spikeTimes = uint64(S0.viTime_spk - 1);
    spikeClusters = uint32(S_clu.viClu);
    channelMap = int32(P.viSite2Chan - 1);
    channelPos = P.mrSiteXY;
    
    % compute 1st three sitewise principal components for each spike, for
    % all neighboring sites
    fprintf('exporting features...');
    pcFeatures = zeros(nSpikes, 3, size(P.miSites, 1), 'single');
    spikeChannels = zeros(nSpikes, size(P.miSites, 1), 'uint32');
    vppFeatures = zeros(nSpikes, 2, size(P.miSites, 1), 'single');

    for iSite = 1:nSites
        siteSpikes = S0.viSite_spk == iSite;
        siteNeighbors = P.miSites(:, iSite);
        spikeWaveforms = tnWav_spk(:, :, siteSpikes);
        
        % PC features
        for jSite = 1:numel(siteNeighbors)
            siteWaveforms = squeeze(spikeWaveforms(:, jSite, :));
            swCentered = siteWaveforms - mean(siteWaveforms, 1);
            covMat = (swCentered * swCentered')/size(swCentered, 2);

            [V, D] = eig(covMat);
            [~, indices] = sort(diag(D), 'descend');
            uT = V(:, indices(1:3));
            pcFeatures(siteSpikes, :, jSite) = siteWaveforms'*uT;
        end

        % vpp (min/max) features
        tnWav_spk1 = jrclust.utils.filtTouV(tnWav_spk_sites_(find(siteSpikes), siteNeighbors, S0), P);
        [mins, maxes] = multifun_(@(x) squeeze(x), min(tnWav_spk1), max(tnWav_spk1));
        vppFeatures(siteSpikes, 1, :) = mins';
        vppFeatures(siteSpikes, 2, :) = maxes';

        spikeChannels(siteSpikes, :) = repmat(P.viSite2Chan(siteNeighbors)' - 1, size(swCentered, 2), 1);
    end
    fprintf('done\n');    

    if ~isempty(savePath)
        writeNPY(spikeTimes, fullfile(savePath, 'spike_times.npy'));
        writeNPY(spikeClusters, fullfile(savePath, 'spike_clusters.npy'));
        writeNPY(channelMap, fullfile(savePath, 'channel_map.npy'));
        writeNPY(channelPos, fullfile(savePath, 'channel_positions.npy'));
        writeNPY(pcFeatures, fullfile(savePath, 'pc_features.npy'));
        writeNPY(vppFeatures, fullfile(savePath, 'vpp_features.npy'));
        writeNPY(spikeChannels, fullfile(savePath, 'spike_channels.npy'));
    end
    
    %make params file
    if ~exist(fullfile(savePath,'params.py'),'file')
        fid = fopen(fullfile(savePath,'params.py'), 'w');
        
        if ~isempty(P.vcFile)
            [~, fname, ext] = fileparts(P.vcFile);
            fprintf(fid, 'dat_path = r''%s%s''\n', fname, ext);
        elseif ischar(P.csFile_merge)
            [~, fname, ext] = fileparts(P.csFile_merge);
            fprintf(fid, 'dat_path = r''%s%s''\n', fname, ext);
        else
            fprintf(fid, 'dat_path = [');
            for i=1:numel(P.csFile_merge) - 1
                fn = P.csFile_merge{i};
                fprintf(fid, 'r''%s'',\n             ', fn);
            end
            fn = P.csFile_merge{end};
            fprintf(fid, 'r''%s'']', fn);
        end
        
        fprintf(fid, 'n_channels_dat = %i\n', P.nChans);
        fprintf(fid, 'n_samples_raw = %i\n', S0.dimm_raw(1));
        fprintf(fid, 'jrc_feature = r''$s''\n', P.vcFet);
        fprintf(fid, 'dtype = r''%s''\n', P.vcDataType);
        fprintf(fid, 'offset = %i\n', P.header_offset);
        fprintf(fid, 'sample_rate = %i\n', P.sRateHz);
        fprintf(fid,'hp_filtered = False');
        fclose(fid);
    end
end

