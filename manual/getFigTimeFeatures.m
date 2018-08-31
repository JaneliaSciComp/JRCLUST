%--------------------------------------------------------------------------
function [features, spikeTimes, yLabel, spikes] = getFigTimeFeatures(site, cluster, S0)
    % just specify site to obtain background info
    % 2016 07 07 JJJ
    % return feature corresponding to a site and cluster
    % requiring subsampled info: cvrVpp_site and cmrFet_site. store in S0

    if nargin < 2
        cluster = [];
    end
    if nargin < 3
        S0 = get(0, 'UserData');
    end

    S_clu = S0.S_clu;
    P = S0.P;
    if ~isfield(P, 'displayFeature')
        P.displayFeature = 'vpp';
    end

    % [features, spikes] = getFet_clu_(site, cluster, S0);
    MAX_SAMPLE = 10000; % max points to display
    spikeSites = S0.spikeSites;

    if isempty(cluster) % select spikes based on sites
        nSitesOfInterest = 1 + round(P.maxSite);
        sitesOfInterest = P.miSites(1:nSitesOfInterest, site);

        try
            spikes = cell2mat_([S0.cviSpk_site(sitesOfInterest)]');
        catch
            spikes = find(ismember(spikeSites, sitesOfInterest));
        end
        spikes = randomSubsample(spikes, MAX_SAMPLE);
    else
        spikes = S_clu.spikesByCluster{cluster};
    end

    switch lower(P.displayFeature)
        case {'vmin', 'vpp'}
            mrWav_spk1 = squeeze_(tnWav2uV_(getSpikeWaveformsSites(spikes, site, S0), P));
            features = max(mrWav_spk1) - min(mrWav_spk1); % >= 0

        case 'kilosort'
            if ~isfield(S0, 'pcTime')
                S0.pcTime = 1;
                set(0, 'UserData', S0);
            end
            
            features = getKilosortFeaturesSites(spikes, site, S0, S0.pcTime);

        case 'cov'
            features = abs(calc_cov_spk_(spikes, site));

        case {'pca', 'gpca'}
            features = abs(pca_pc_spk_(spikes, site));

        case {'ppca', 'private pca'}
            [mrPv1, mrPv2] = pca_pv_clu_(site, S0.primarySelectedCluster);
            features = abs(pca_pc_spk_(spikes, site, mrPv1, mrPv2));

        otherwise
            error('not implemented yet');
    end

    features = squeeze(features);
    % end getFet_clu_
    spikeTimes = double(S0.spikeTimes(spikes)) / P.sampleRateHz;

    % label
    switch lower(P.displayFeature)
        case {'vpp', 'vmin'} %voltage feature
            yLabel = sprintf('Site %d (\\mu%s)', site, P.displayFeature);
            
        case 'kilosort'
            yLabel = sprintf('Site %d PC%d', site, S0.pcTime);

        otherwise %other feature options
            yLabel = sprintf('Site %d (%s)', site, P.displayFeature);
    end

end % function
