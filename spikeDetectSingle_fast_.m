%--------------------------------------------------------------------------
% 11/15/17 JJJ: Cast the threshold like the vrWav1
function [viSpk1, vrSpk1, thresh1] = spikeDetectSingle_fast_(vrWav1, P, thresh1)
    % P: spkThresh, qqSample, qqFactor, useGPU, uV_per_bit
    % vrWav1 can be either single or int16
    % 6/27/17 JJJ: bugfix: hard set threshold is applied

    % Determine threshold
    MAX_SAMPLE_QQ = 300000;
    % fSpikeRefrac_site = 0;
    if nargin < 3, thresh1 = []; end
    if nargin < 2, P = struct('spkThresh', [], 'qqFactor', 5); end
    if ~isempty(get_(P, 'spkThresh')), thresh1 = P.spkThresh; end

    if thresh1==0, [viSpk1, vrSpk1] = deal([]); return; end % bad site
    if isempty(thresh1)
        thresh1 = median(abs(subsample_vr_(vrWav1, MAX_SAMPLE_QQ)));
        thresh1 = single(thresh1)* P.qqFactor / 0.6745;
    end
    thresh1 = cast(thresh1, 'like', vrWav1); % JJJ 11/5/17

    % detect valley turning point. cannot detect bipolar
    % pick spikes crossing at least three samples
    nneigh_min = get_set_(P, 'nneigh_min_detect', 0);
    viSpk1 = find_peak_(vrWav1, thresh1, nneigh_min);
    if get_set_(P, 'fDetectBipolar', 0)
        viSpk1 = [viSpk1; find_peak_(-vrWav1, thresh1, nneigh_min)];
        viSpk1 = sort(viSpk1);
    end
    if isempty(viSpk1)
        viSpk1 = double([]);
        vrSpk1 = int16([]);
    else
        vrSpk1 = vrWav1(viSpk1);
        % Remove spikes too large
        spkThresh_max_uV = get_set_(P, 'spkThresh_max_uV', []);
        if ~isempty(spkThresh_max_uV)
            thresh_max1 = abs(spkThresh_max_uV) / get_set_(P, 'uV_per_bit', 1);
            thresh_max1 = cast(thresh_max1, 'like', vrSpk1);
            viA1 = find(abs(vrSpk1) < abs(thresh_max1));
            viSpk1 = viSpk1(viA1);
            vrSpk1 = vrSpk1(viA1);
        end
    end

    % apply spike merging on the same site
    % nRefrac = int32(abs(P.spkRefrac));
    % if P.refrac_factor > 1
    %     nRefrac = int32(round(double(nRefrac) * P.refrac_factor));
    % end
    if isGpu_(viSpk1)
        [viSpk1, vrSpk1, thresh1] = multifun_(@gather, viSpk1, vrSpk1, thresh1);
    end
    % if fSpikeRefrac_site %perform spike refractive period per site (affects exact mode)
    %     [viSpk1, vrSpk1] = spike_refrac_(viSpk1, vrSpk1, [], nRefrac); %same site spikes
    % end
end %func
