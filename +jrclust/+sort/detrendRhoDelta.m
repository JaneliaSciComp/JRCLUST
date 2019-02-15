function [centers, logRho, zscores] = detrendRhoDelta(hClust, spikesBySite, fLocal, hCfg)
    %DETREND Detrend rho-delta plot to identify cluster centers (high rho, high delta)
    logRho = log10(hClust.spikeRho);
    delta = hClust.spikeDelta;

    if fLocal % detrend for each site
        rhosBySite = cellfun(@(spikes) hClust.spikeRho(spikes), spikesBySite, 'UniformOutput', 0);
        deltasBySite = cellfun(@(spikes) delta(spikes), spikesBySite, 'UniformOutput', 0);

        centersBySite = cell(size(spikesBySite));
        zscores = zeros(size(delta), 'like', delta);

        for iSite = 1:numel(spikesBySite)
            siteSpikes = spikesBySite{iSite};
            if isempty(siteSpikes)
                continue;
            end

            logRhoSite = log10(rhosBySite{iSite});
            deltaSite = deltasBySite{iSite};

            % select only those rho values between 10^cutoff and 0.1
            % and delta values between 0 and 1
            detIndices = find(logRhoSite > hCfg.log10RhoCut & logRhoSite < -1 & deltaSite > 0 & deltaSite < 1 & isfinite(logRhoSite) & isfinite(deltaSite));
            [yhat, zsite] = detrendQuadratic(logRhoSite, deltaSite, detIndices);

            % spikes with exact same features
            yhat(deltaSite == 0) = nan;
            zsite(deltaSite == 0) = nan;

            % get the indices of the maxClustersSite largest detrended
            % values where the rho values beat the cutoff
            centers_ = nLargest(yhat, hCfg.maxClustersSite, find(logRhoSite > hCfg.log10RhoCut & ~isnan(yhat)));
            if isempty(centers_)
                continue;
            end

            centersBySite{iSite} = siteSpikes(centers_);
            zscores(siteSpikes) = zsite;
        end

        centers = cell2vec(centersBySite);
    else % detrend globally
        detIndices = find(delta > 0 & delta < 1 & hClust.spikeRho > 10^hCfg.log10RhoCut & hClust.spikeRho < .1 & isfinite(logRho) & isfinite(delta));
        [~, zscores] = detrendQuadratic(logRho, delta, detIndices);

        zscores(delta == 0) = nan;
        [centers, zLargest] = nLargest(zscores, hCfg.maxClustersSite*numel(spikesBySite), find(hClust.spikeRho > 10^hCfg.log10RhoCut & ~isnan(zscores)));
        centers(zLargest < 10^hCfg.log10DeltaCut) = [];
    end
end

%% LOCAL FUNCTIONS
function [yhat, zscores] = detrendQuadratic(x, y, indices)
    %DETRENDQUADRATIC Remove quadratic component from y~x
    x = x(:);
    y = y(:);

    eps_ = eps(class(x));
    X = [x.^2 + eps_, x + eps_, ones(numel(x), 1)];
    b = X(indices, :)\y(indices);

    yhat = y - X*b;

    zscores = (yhat - nanmean(yhat(indices))) / nanstd(yhat(indices));
end

function  [iLargest, largest] = nLargest(vals, n, indices)
    %NLARGEST Get the (indices of the) n largest elements in an array
    n = min(n, numel(indices));
    [~, argsort] = sort(vals(indices), 'descend');

    iLargest = indices(argsort(1:n));
    largest = vals(iLargest);
end

function vals = cell2vec(vals)
    %CELL2VEC Convert cell array to vector
    vals = vals(:);
    vals = vals(cellfun(@(x) ~isempty(x), vals));

    for i = 1:numel(vals)
        v = vals{i};
        vals{i} = v(:);
    end

    vals = cell2mat(vals);
end
