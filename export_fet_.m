%--------------------------------------------------------------------------
function export_fet_(P)
    % export feature matrix to workspace
    S0 = load(strrep(P.paramFile, '.prm', '_jrc.mat'));
    trFet = load_bin_(strrep(P.paramFile, '.prm', '_features.bin'), 'single', S0.featureDims);
    assignWorkspace_(trFet);
end %func
