%--------------------------------------------------------------------------
function export_spkwav_(P, vcArg2, fDiff)
    % Export spike waveforms organized by clusters
    % export_spkwav_(P)
    % export_spkwav_(P, viClu)
    global fDebug_ui

    if nargin<2, vcArg2 = ''; end
    if nargin<3, fDiff = 0; end
    if isempty(vcArg2)
        viClu = [];
        fPlot = 0;
    else
        viClu = str2num(vcArg2);
        if isempty(viClu), fprintf(2, 'Invalid Cluster #: %s', vcArg2); end
        fPlot = numel(viClu)==1;
    end

    % Load data
    S0 = load_cached_(P);
    if isempty(S0), fprintf(2, 'Not clustered yet'); return; end
    if ~isfield(S0, 'S_clu'), fprintf(2, 'Not clustered yet'); return; end
    S_clu = S0.S_clu;

    % Collect waveforms by clusters
    ctrWav_clu = cell(1, S_clu.nClu);
    miSite_clu = P.miSites(:,S_clu.viSite_clu);
    fprintf('Collecting spikes from clusters\n\t'); t1=tic;
    if isempty(viClu), viClu = 1:S_clu.nClu; end
    for iClu = viClu
        tnWav_clu1 = tnWav_spk_sites_(S_clu.cviSpk_clu{iClu}, miSite_clu(:,iClu), S0);
        if fDiff
            ctrWav_clu{iClu} = tnWav_clu1;
        else
            ctrWav_clu{iClu} = tnWav2uV_(tnWav_clu1, P);
        end
        %     ctrWav_clu{iClu} = (cumsum(bit2uV_(meanSubt_(tnWav_clu1))));
        %     ctrWav_clu{iClu} = (bit2uV_(meanSubt_(tnWav_clu1)));
        fprintf('.');
    end
    fprintf('\n\ttook %0.1fs\n', toc(t1));

    if fPlot
        iClu = viClu;
        nT_spk = (diff(P.spkLim)+1);
        nSpk1 = S_clu.vnSpk_clu(iClu);
        hFig = createFigure(sprintf('Fig_clu%d', iClu), [0 0 .5 1], P.prmFile, 1, 1);
        multiplot([], P.maxAmp, [], ctrWav_clu{iClu}, miSite_clu(:,iClu));
        xlabel('Spike #'); ylabel('Site #');
        set(gca, 'YTick', range_(miSite_clu(:,iClu)), ...
        'XTick', (1:nSpk1) * nT_spk - P.spkLim(2), 'XTickLabel', 1:nSpk1);
        axis_([0, nSpk1*nT_spk+1, (limit_(miSite_clu(:,iClu)) + [-1,1])]);
        grid on;
        title(sprintf('Cluster%d (n=%d), Scale=%0.1f uV', iClu, nSpk1, P.maxAmp));
        mouse_figure(hFig);
        eval(sprintf('trWav_clu%d = ctrWav_clu{iClu};', iClu));
        eval(sprintf('viSites_clu%d = miSite_clu(:,iClu);', iClu));
        eval(sprintf('assignWorkspace_(trWav_clu%d, viSites_clu%d);', iClu, iClu));
        %     if get_set_([], 'fDebug_ui', 0), close_(hFig); end
        if fDebug_ui==1, close_(hFig); end
    else
        assignWorkspace_(ctrWav_clu, miSite_clu);
    end
end %func
