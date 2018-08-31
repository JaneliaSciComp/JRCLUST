%--------------------------------------------------------------------------
function [featuresMin, featuresMax] = getFeatureForSpikes(spikes, sitesOfInterest, S0)
    % get feature for the spikes of interest (TODO: this is duplicated
    % effort; delete)

    if nargin < 3
        S0 = get(0, 'UserData');
    end

    P = S0.P;

    switch lower(P.displayFeature)
        case {'vmin', 'vpp'}
            spikeWaveforms = getSpikeWaveformsSites(spikes, sitesOfInterest, S0);
            spikeWaveforms = tnWav2uV_(spikeWaveforms, P);

            featuresMin = squeeze(abs(min(spikeWaveforms))); % nSites x nSpikes
            featuresMax = squeeze(abs(max(spikeWaveforms)));

        case {'cov', 'spacetime'}
            [featuresMin, featuresMax] = calc_cov_spk_(spikes, sitesOfInterest);

        case 'pca'
            [featuresMin, featuresMax] = pca_pc_spk_(spikes, sitesOfInterest);

        case 'kilosort'
            [featuresMin, featuresMax] = getKilosortFeaturesSites(spikes, sitesOfInterest, S0, S0.pcPair);

        otherwise
            error('not implemented yet');
    end
end % function
