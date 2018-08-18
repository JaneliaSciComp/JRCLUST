%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function spikeFeatures = getSpikeFeatures(P)
    if nargin < 1
        P = get0_('P');
    end

    spikeFeatures = load_bin_(strrep(P.paramFile, '.prm', '_features.bin'), 'single', get0_('featureDims'));
end %func
