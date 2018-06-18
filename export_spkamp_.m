%--------------------------------------------------------------------------
function export_spkamp_(P, vcArg2)
    % export spike amplitudes (Vpp, Vmin, Vmax) in uV to workspace, organize by clusters

    if nargin<2, vcArg2 = ''; end
    if isempty(vcArg2)
        viClu = [];
        fSingleUnit = 0;
    else
        viClu = str2num(vcArg2);
        if isempty(viClu), fprintf(2, 'Invalid Cluster #: %s', vcArg2); end
        fSingleUnit = numel(viClu)==1;
    end

    % Load data
    S0 = load_cached_(P);
    if isempty(S0), fprintf(2, 'Not clustered yet'); return; end
    if ~isfield(S0, 'S_clu'), fprintf(2, 'Not clustered yet'); return; end
    S_clu = S0.S_clu;

    % Collect waveforms by clusters
    [cmrVpp_clu, cmrVmin_clu, cmrVmax_clu] = deal(cell(1, S_clu.nClu));
    miSite_clu = P.miSites(:,S_clu.viSite_clu);
    fprintf('Calculating spike amplitudes from clusters\n\t'); t1=tic;
    if isempty(viClu), viClu = 1:S_clu.nClu; end
    for iClu = viClu
        trWav_clu1 = tnWav2uV_(tnWav_spk_sites_(S_clu.cviSpk_clu{iClu}, miSite_clu(:,iClu), S0), P);
        cmrVmin_clu{iClu} = shiftdim(min(trWav_clu1));
        cmrVmax_clu{iClu} = shiftdim(max(trWav_clu1));
        cmrVpp_clu{iClu} = cmrVmax_clu{iClu} - cmrVmin_clu{iClu};
        fprintf('.');
    end
    fprintf('\n\ttook %0.1fs\n', toc(t1));

    % export to workspace
    if fSingleUnit
        iClu = viClu;
        eval(sprintf('mrVmax_clu%d = cmrVmax_clu{iClu};', iClu));
        eval(sprintf('mrVmin_clu%d = cmrVmin_clu{iClu};', iClu));
        eval(sprintf('mrVpp_clu%d = cmrVpp_clu{iClu};', iClu));
        eval(sprintf('viSites_clu%d = miSite_clu(:,iClu);', iClu));
        eval(sprintf('assignWorkspace_(mrVmax_clu%d, mrVmin_clu%d, mrVpp_clu%d, viSites_clu%d);', iClu, iClu, iClu, iClu));
    else
        assignWorkspace_(cmrVpp_clu, cmrVmax_clu, cmrVmin_clu, miSite_clu);
    end
end %func
