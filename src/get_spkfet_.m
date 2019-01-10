%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function trWav_fet = get_spkfet_(P)
    trWav_fet = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkfet.jrc'), 'single', get0_('dimm_fet'));
end %func
