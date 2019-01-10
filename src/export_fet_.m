%--------------------------------------------------------------------------
function export_fet_(P)
    % export feature matrix to workspace
    S0 = load(jrclust.utils.subsExt(P.vcFile_prm, '_jrc.mat'));
    trFet = jrclust.utils.readBin(jrclust.utils.subsExt(P.vcFile_prm, '_spkfet.jrc'), 'single', S0.dimm_fet);
    assignWorkspace_(trFet);
end %func
