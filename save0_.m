%--------------------------------------------------------------------------
% 11/10/17: Removed recursive saving
% 9/29/17 Updating the version number when saving JRCLUST
function S0 = save0_(vcFile_mat, fSkip_fig)
    if nargin<2, fSkip_fig = 0; end
    % save S0 structure to a mat file
    try
        fprintf('Saving UserData to %s...\n', vcFile_mat);
        warning off;
        S0 = get(0, 'UserData'); %add gather script
        if isfield(S0, 'S0'), S0 = rmfield(S0, 'S0'); end % Remove recursive saving

        % update version number
        S0.P.version = jrcVersion();
        P = S0.P;
        setUserData(P);

        struct_save_(S0, vcFile_mat, 1);
        vcFile_prm = S0.P.paramFile;
        exportParams(vcFile_prm, strrep(vcFile_prm, '.prm', '_full.prm'), 0);

        % save the rho-delta plot
        if fSkip_fig, return; end
        if ~isfield(S0, 'S_clu') || ~get_set_(P, 'fSavePlot_RD', 1), return; end
        try
            if ~isfield(S0, 'clusterCenters'), return; end % skip kilosort
            save_fig_(strrep(P.paramFile, '.prm', '_RD.png'), plot_rd_(P, S0), 1);
            fprintf('\tYou can use ''jrc plot-rd'' command to plot this figure.\n');
        catch
            fprintf(2, 'Failed to save the rho-delta plot: %s.\n', lasterr());
        end
    catch
        disperr_();
    end
end %func
