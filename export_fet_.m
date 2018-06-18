%--------------------------------------------------------------------------
function export_fet_(P)
    % export feature matrix to workspace
    S0 = load(strrep(P.vcFile_prm, '.prm', '_jrc.mat'));
    trFet = load_bin_(strrep(P.vcFile_prm, '.prm', '_spkfet.jrc'), 'single', S0.dimm_fet);
    assignWorkspace_(trFet);
end %func
