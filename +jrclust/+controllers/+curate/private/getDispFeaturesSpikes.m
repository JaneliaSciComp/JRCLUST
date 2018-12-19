%--------------------------------------------------------------------------
%function [mrMin, mrMax] = getFet_spk_(spikes, sites, S0)
function [mrMin, mrMax] = getDispFeaturesSpikes(hClust, spikes, sites)
    %GETDISPFEATURESSPIKES Get display feature for the spikes of interest
    hCfg = hClust.hCfg;

    switch lower(hCfg.vcFet_show)
        case {'vmin', 'vpp'}
            tnWav_spk1 = jrclust.utils.filtTouV(jrclust.utils.getSampledWindows(hClust, spikes, sites), hCfg);
            [mrMin, mrMax] = multifun_(@(x)abs(permute(x,[2,3,1])), min(tnWav_spk1), max(tnWav_spk1));

        case {'cov', 'spacetime'}
            [mrMin, mrMax] = getSpikeCov(hClust, spikes, sites);

        case 'pca'
            [mrMin, mrMax] = pca_pc_spk_(spikes, sites); % get all spikes whose center lies in certain range

        case 'kilosort'
            [mrMin, mrMax] = ks_fet_spk_(spikes, sites, S0);

        otherwise
            error('not implemented yet');
    end
end %func
