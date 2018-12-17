function centers = regressCenters(sRes, spikesBySite, thresh)
    % Calculate regression line to select delta and rho
    % Hidehiko Inagako, 20160425

    if nargin < 2
        thresh = 0.2;
    end

    fDebug = 0;
    fDebug1 = 0;

    log10Rho = log10(sRes.spikeRho);
    log10Delta = log10(sRes.spikeDelta);

    xRange = -4:0.1:1.5;
    %xRange = floor(min(log10Rho)):0.1:ceil(max(log10Rho));
    yRange = -3.5:0.1:1.5;
    %yRange = floor(min(log10Delta)):0.1:ceil(max(log10Delta));

    nSites = numel(spikesBySite);

    siteCenters = cell(nSites, 1);

    for iSite = 1:nSites    
        siteSpikes = spikesBySite{iSite};
        if isempty(siteSpikes)
            continue;
        end

        siteDelta = log10Delta(siteSpikes);
        siteRho = log10Rho(siteSpikes);

        % acquire density of scatter 
        I2 = hist3([siteRho(:), siteDelta(:)], 'edges', {xRange, yRange});
        I2 = medfilt2(I2');

        XData = xRange;
        YData = nan(size(xRange));

        for iX = 1:size(I2, 2) % find largest delta bin for this rho bin
            idxFirst = find(I2(:, iX) > 1, 1, 'last');
            if ~isempty(idxFirst)
                YData(iX) = yRange(idxFirst);
            end
        end

        emptyCols = isnan(YData);
        XData(emptyCols) = [];
        YData(emptyCols) = [];

        % regression of upper edge line
        X = [ones(length(XData), 1) XData'];
        b = X \ YData';
        vrLDF1 = b(1) + b(2)*siteRho + thresh; % fit line
        siteCenters{iSite} = siteSpikes(siteDelta > vrLDF1);

        if fDebug1
            figure; hold on;
            imagesc(I2, 'XData', xRange, 'YData', yRange);
            axis xy;
            plot(siteRho, siteDelta, '.');
            scent = ismember(siteSpikes, siteCenters{iSite});
            plot(siteRho(scent), siteDelta(scent), '*')
            plot(XData, YData, 'r-');
            plot(siteRho, vrLDF1, 'k-');
        end
    end

    centers = jrclust.utils.neCell2mat(siteCenters);

    if fDebug
        figure; plot(log10Rho(1:10:end), log10Delta(1:10:end), '.', log10Rho(centers), log10Delta(centers), 'r.');
    end
end

