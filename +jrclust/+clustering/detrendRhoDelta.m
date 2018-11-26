function [centers, x, zscores] = detrendRhoDelta(clusterData, spikesBySite, fLocal, hCfg)
    %DETREND Detrend rho-delta plot to identify cluster centers (high rho, high delta)

    x = log10(clusterData.spikeRho);

    % detrend for each site and apply
    if fLocal
        rhosBySite = cellfun(@(spikes) clusterData.spikeRho(spikes), spikesBySite, 'UniformOutput', 0);
        deltasBySite = cellfun(@(spikes) clusterData.spikeDelta(spikes), spikesBySite, 'UniformOutput', 0);

        centersBySite = cell(size(spikesBySite));
        zscores = zeros(size(clusterData.spikeDelta), 'like', clusterData.spikeDelta);

        for iSite = 1:numel(spikesBySite)
            siteSpikes = spikesBySite{iSite};
            if isempty(siteSpikes)
                continue;
            end

            rhoSite = rhosBySite{iSite};
            deltaSite = deltasBySite{iSite};

            log10RhoSite = log10(rhoSite);

            % select only those rho values between 10^cutoff and 0.1
            % and delta values between 0 and 1
            detIndices = find(log10RhoSite > hCfg.log10RhoCut & log10RhoSite < -1 & deltaSite > 0 & deltaSite < 1 & isfinite(log10RhoSite) & isfinite(deltaSite));
            [yhat, zsite] = detrendQuadratic(log10RhoSite(detIndices), deltaSite(detIndices));

            yhat(deltaSite == 0) = nan;
            zsite(deltaSite == 0) = nan;

            % get the indices of the maxClustersSite largest detrended
            % values where the rho values beat the cutoff
            centers_ = nLargest(yhat, hCfg.maxClustersSite, find(log10RhoSite > hCfg.log10RhoCut & ~isnan(yhat)));
            if isempty(centers_)
                continue;
            end

            centersBySite{iSite} = siteSpikes(centers_);
            zscores(siteSpikes) = zsite;
        end

        centers = cell2vec(centersBySite);
    else
        y = clusterData.spikeDelta;
        detIndices = find(clusterData.spikeDelta < 1 & clusterData.spikeDelta > 0 & clusterData.spikeRho > 10^hCfg.rho_cut & clusterData.spikeRho < .1 & isfinite(x) & isfinite(y));
        [~, zscores] = detrendQuadratic(x(detIndices), y(detIndices));

        zscores(clusterData.spikeDelta == 0) = nan;
        [centers, zLargest] = nLargest(zscores, hCfg.maxClustersSite*numel(spikesBySite), find(clusterData.spikeRho > 10^hCfg.log10RhoCut & ~isnan(zscores)));
        centers(zLargest < 10^hCfg.log10DeltaCut) = [];
    end
end

%% LOCAL FUNCTIONS
function [yhat, zscores] = detrendQuadratic(x, y)
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
