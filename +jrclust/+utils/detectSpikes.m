function [viTime_spk, vrAmp_spk, viSite_spk] = detectSpikes(samplesIn, siteThresh, keepMe, hCfg)
    %DETECTSPIKES
    if isfield(hCfg, 'fMerge_spk')
        fMerge_spk = hCfg.fMerge_spk;
    else
        fMerge_spk = true;
    end

    [n1, nSites, ~] = size(samplesIn);
    [cviSpk_site, cvrSpk_site] = deal(cell(nSites, 1));

    fprintf('\tDetecting spikes from each channel.\n\t\t');
    td = tic;
    for iSite = 1:nSites
        % Find spikes
        [viSpk11, vrSpk11] = spikeDetectSingle_fast_(samplesIn(:, iSite), hCfg, siteThresh(iSite));
        fprintf('.');

        % Reject global mean
        if isempty(keepMe)
            cviSpk_site{iSite} = viSpk11;
            cvrSpk_site{iSite} = vrSpk11;
        else
            [cviSpk_site{iSite}, cvrSpk_site{iSite}] = select_vr_(viSpk11, vrSpk11, find(keepMe(viSpk11)));
        end
    end

    siteThresh = jrclust.utils.tryGather(siteThresh);
    nSpks1 = sum(cellfun(@numel, cviSpk_site));
    fprintf('\n\t\tDetected %d spikes from %d sites; took %0.1fs.\n', nSpks1, nSites, toc(td));

    % Group spiking events using vrWav_mean1. already sorted by time
    if fMerge_spk
        fprintf('\tMerging spikes...');
        t2=tic;
        [viTime_spk, vrAmp_spk, viSite_spk] = spikeMerge_(cviSpk_site, cvrSpk_site, hCfg);
        fprintf('\t%d spiking events found; took %0.1fs\n', numel(viSite_spk), toc(t2));
    else
        viTime_spk = jrclust.utils.neCell2mat(cviSpk_site);
        vrAmp_spk = jrclust.utils.neCell2mat(cvrSpk_site);
        viSite_spk = cell2vi_(cviSpk_site);

        % sort by time
        [viTime_spk, viSrt] = sort(viTime_spk, 'ascend');
        [vrAmp_spk, viSite_spk] = multifun_(@(x)x(viSrt), vrAmp_spk, viSite_spk);
    end
    vrAmp_spk = jrclust.utils.tryGather(vrAmp_spk);

    % Group all sites in the same shank
    if hCfg.fGroup_shank'
        [viSite_spk] = group_shank_(viSite_spk, hCfg); % change the site location to the shank center
    end
end

%% LOCAL FUNCTIONS
function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_fast_(vrWav1, hCfg, thresh1)
    % hCfg: spkThresh, qqSample, qqFactor, fGpu, uV_per_bit
    % vrWav1 can be either single or int16
    % 6/27/17 JJJ: bugfix: hard set threshold is applied

    % Determine threshold
    MAX_SAMPLE_QQ = 300000;
    % fSpikeRefrac_site = 0;
    if nargin < 3
        thresh1 = [];
    end
    if nargin < 2
        hCfg = struct('spkThresh', [], 'qqFactor', 5);
    end

    if ~isempty(get_(hCfg, 'spkThresh'))
        thresh1 = hCfg.spkThresh;
    end

    if thresh1 == 0
        [viSpk1, vrSpk1] = deal([]);
        return;
    end % bad site

    if isempty(thresh1)
        thresh1 = median(abs(subsample_vr_(vrWav1, MAX_SAMPLE_QQ)));
        thresh1 = single(thresh1)* hCfg.qqFactor / 0.6745;
    end
    thresh1 = cast(thresh1, 'like', vrWav1); % JJJ 11/5/17

    % detect valley turning point. cannot detect bipolar
    % pick spikes crossing at least three samples
    viSpk1 = find_peak_(vrWav1, thresh1, hCfg.nneigh_min_detect);
    if hCfg.fDetectBipolar
        viSpk1 = [viSpk1; find_peak_(-vrWav1, thresh1, hCfg.nneigh_min_detect)];
        viSpk1 = sort(viSpk1);
    end
    if isempty(viSpk1)
        viSpk1 = double([]);
        vrSpk1 = int16([]);
    else
        vrSpk1 = vrWav1(viSpk1);
        % Remove spikes too large
        if ~isempty(hCfg.spkThresh_max_uV)
            thresh_max1 = abs(hCfg.spkThresh_max_uV) / hCfg.bitScaling;
            thresh_max1 = cast(thresh_max1, 'like', vrSpk1);
            viA1 = find(abs(vrSpk1) < abs(thresh_max1));
            viSpk1 = viSpk1(viA1);
            vrSpk1 = vrSpk1(viA1);
        end
    end

    % apply spike merging on the same site
    % nRefrac = int32(abs(hCfg.spkRefrac));
    % if hCfg.refrac_factor > 1
    %     nRefrac = int32(round(double(nRefrac) * hCfg.refrac_factor));
    % end
    if isa(viSpk1, 'gpuArray')
        [viSpk1, vrSpk1, thresh1] = multifun_(@gather, viSpk1, vrSpk1, thresh1);
    end
end

function varargout = select_vr_(varargin)
    % [var1, var2, ...] = select_vr(var1, var2, ..., index)

    % sort ascend
    viKeep = varargin{end};
    if islogical(viKeep), viKeep = find(viKeep); end
    for i=1:(nargin-1)
        if isvector(varargin{i})
            varargout{i} = varargin{i}(viKeep);
        else
            varargout{i} = varargin{i}(viKeep, :);
        end
    end
end

function [viSpk, vrSpk, viSite] = spikeMerge_(cviSpk, cvrSpk, hCfg)
    % provide spike index (cviSpk) and amplitudes (cvrSPk) per sites

    nSites = numel(cviSpk);
    viSpk = jrclust.utils.neCell2mat(cviSpk);      vrSpk = jrclust.utils.neCell2mat(cvrSpk);
    viSite = jrclust.utils.neCell2mat(cellfun(@(vi,i)repmat(i,size(vi)), cviSpk, num2cell((1:nSites)'), 'UniformOutput', false));
    [viSpk, viSrt] = sort(viSpk);   vrSpk = vrSpk(viSrt);   viSite = viSite(viSrt);
    viSite = int32(viSite);
    viSpk = int32(viSpk);

    [cviSpkA, cvrSpkA, cviSiteA] = deal(cell(nSites,1));

    try
        parfor iSite = 1:nSites %parfor speedup: 2x %parfor
            try
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, hCfg);
            catch
                disperr_();
            end
        end
    catch
        for iSite = 1:nSites
            try
                [cviSpkA{iSite}, cvrSpkA{iSite}, cviSiteA{iSite}] = ...
                spikeMerge_single_(viSpk, vrSpk, viSite, iSite, hCfg);
            catch
                disperr_();
            end
        end
    end

    % merge parfor output and sort
    viSpk = jrclust.utils.neCell2mat(cviSpkA);
    vrSpk = jrclust.utils.neCell2mat(cvrSpkA);
    viSite = jrclust.utils.neCell2mat(cviSiteA);
    [viSpk, viSrt] = sort(viSpk); %sort by time
    vrSpk = jrclust.utils.tryGather(vrSpk(viSrt));
    viSite = viSite(viSrt);
end

function [viSite_spk] = group_shank_(viSite_spk, hCfg)
    nSites = numel(hCfg.siteMap);
    site2site = zeros([nSites, 1], 'like', viSite_spk);
    [~, b, c] = unique(hCfg.viShank_site);
    site2site(hCfg.siteMap) = b(c);
    viSite_spk = site2site(viSite_spk);
end

