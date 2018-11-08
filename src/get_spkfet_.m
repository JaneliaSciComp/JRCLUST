%--------------------------------------------------------------------------
% 10/11/17 JJJ: created
function trWav_fet = get_spkfet_(P)
    trWav_fet = jrclust.utils.readBin(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', get0_('dimm_fet'));
end %func
