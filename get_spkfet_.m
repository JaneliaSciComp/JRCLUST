%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function trWav_fet = get_spkfet_(P)
    trWav_fet = load_bin_(strrep(P.paramFile, '.prm', '_features.bin'), 'single', get0_('featureDims'));
end %func
