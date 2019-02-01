function varargout = spikeCov(spikeWindows, hCfg)
    %SPIKECOV Compute covariance feature for spikes
    %   spikeWindows: nSamples x nSpikes x nSites
    vnDelay_fet = hCfg.getOr('vnDelay_fet', [0, 3]);

    [nSamples, ~, nSites] = size(spikeWindows);
    [shiftsBack, shiftsForward] = jrclust.utils.shiftRange(nSamples, [], vnDelay_fet);
    cmrFet = cell(numel(vnDelay_fet), 1);

    for iDelay = 1:numel(vnDelay_fet)
        mr1_ = jrclust.utils.meanSubtract(spikeWindows(shiftsForward{iDelay}, :, 1));
        mr1_ = bsxfun(@rdivide, mr1_, sqrt(mean(mr1_.^2))); %zscore fast

        tr1_ = repmat(mr1_, [1,1,nSites]);
        tr2_ = jrclust.utils.meanSubtract(spikeWindows(shiftsBack{iDelay},:,:));
        cmrFet{iDelay} = permute(mean(tr1_ .* tr2_, 1), [3,2,1]);
    end %for

    if nargout==1
        varargout{1} = jrclust.utils.neCell2mat(cmrFet);
    else
        varargout = cmrFet;
    end
end
