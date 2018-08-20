%--------------------------------------------------------------------------
function [vrY, vrX] = tr2plot_(spikeWaveforms, clusterID, sitesToDisplay, maxAmp, P)
    if nargin < 2
        clusterID = 1;
    end
    
    clusterID = double(clusterID);
    if nargin < 5
        P = get0_('P');
    end

    if nargin < 3
        sitesToDisplay = [];
    end

    if isempty(sitesToDisplay)
        sitesToDisplay = 1:size(spikeWaveforms,2);
    end

    [nSamples, ~, nSpikes] = size(spikeWaveforms);
    nSitesToDisplay = numel(sitesToDisplay);
    spikeWaveforms = single(spikeWaveforms) / maxAmp;
    
    % offset each waveform by its site index
    spikeWaveforms = spikeWaveforms + repmat(single(sitesToDisplay(:)'), [nSamples, 1, nSpikes]);
    spikeWaveforms([1, end], :, :) = nan; % why? -- acl
    vrY = spikeWaveforms(:);

    if nargout >= 2
        vrX = wav_clu_x_(clusterID, P);
        vrX = repmat(vrX(:), [1, nSitesToDisplay * nSpikes]);
        vrX = single(vrX(:));
    end
end
