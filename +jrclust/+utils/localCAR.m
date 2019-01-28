function [spikeWindows, carTraces] = localCAR(spikeWindows, hCfg, nSitesEvt, refSites)
    %LOCALCAR Subtract local common-average reference from spike windows
    %   spikeWindows: nSamples x nSpikes x nSites, single
    if nargin < 3
        nSitesEvt = [];
    end
    if nargin < 4 % usually (always?) the case
        refSites = [];
    end

    if isempty(nSitesEvt)
        nSitesEvt = hCfg.nSitesEvt;
    end

    if nSitesEvt == 1
        carTraces = [];
        return;
    end

    if strcmp(hCfg.getOr('vcSpkRef', 'nmean'), 'nmean') % use n sites having the least SD as reference sites
        if isempty(refSites)
            farHalf = ceil(size(spikeWindows, 3)/2):size(spikeWindows, 3);
            carTraces = mean(spikeWindows(:, :, farHalf), 3);
        else
            innerTraces = spikeWindows(:, :, 1:nSitesEvt);
            innerTraces(:, :, 1) = 0;

            for iSpk1 = 1:numel(refSites)
                innerTraces(:, iSpk1, refSites(iSpk1)) = 0;
            end

            carTraces = sum(innerTraces, 3) / (nSitesEvt-2);
        end
    else
        carTraces = [];
    end

    spikeWindows = spikeWindows(:, :, 1:nSitesEvt);

    shape = size(spikeWindows);
    if ismatrix(spikeWindows)
        shape(end+1) = 1;
    end

    if ~isempty(carTraces) % subtract carTraces from spikeWindows
        spikeWindows = jrclust.utils.meanSubtract(reshape(bsxfun(@minus, reshape(spikeWindows, [], shape(3)), carTraces(:)), shape));
    else
        spikeWindows = jrclust.utils.meanSubtract(spikeWindows);
    end
end
