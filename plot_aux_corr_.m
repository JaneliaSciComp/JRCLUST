%--------------------------------------------------------------------------
% 9/19/17 JJJ: Created for SPARC
function plot_aux_corr_(mrRate_clu, vrWav_aux, vrCorr_aux_clu, vrTime_aux, iCluPlot)
    if nargin<5, iCluPlot = []; end
    % show the firing rate and plot the
    [vrCorr_srt, viSrt] = sort(vrCorr_aux_clu, 'descend');
    nClu = numel(vrCorr_aux_clu);
    [P, S_clu] = get0_('P', 'S_clu');
    P = loadParams(P.paramFile);
    nClu_show = min(get_set_(P, 'nClu_show_aux', 4), nClu);
    vcLabel_aux = get_set_(P, 'vcLabel_aux', 'aux');
    nSubsample_aux = get_set_(P, 'nSubsample_aux', 100);
    if ~isempty(iCluPlot)
        nClu_show = 1;
        viSrt = iCluPlot;
    end

    hFig = createFigure('FigAux', [.5 0 .5 1], P.paramFile,1,1);
    hTabGroup = uitabgroup(hFig);
    for iClu1 = 1:nClu_show
        iClu = viSrt(iClu1);
        htab1 = uitab(hTabGroup, 'Title', sprintf('Clu %d', iClu), 'BackgroundColor', 'w');
        ax_ = axes('Parent', htab1);
        subplot(2, 1, 1);
        ax_ = plotyy(vrTime_aux, mrRate_clu(:,iClu), vrTime_aux, vrWav_aux);
        xlabel('Time (s)');
        ylabel(ax_(1),'Firing Rate (Hz)');
        ylabel(ax_(2), vcLabel_aux);
        iSite_ = S_clu.clusterSites(iClu);
        vcTitle_ = sprintf('Clu %d (Site %d, Chan %d): Corr=%0.3f', ...
        iClu, iSite_, P.chanMap(iSite_), vrCorr_aux_clu(iClu));
        title(vcTitle_);
        set(ax_, 'XLim', vrTime_aux([1,end]));
        grid on;

        subplot(2, 1, 2);
        plot(vrWav_aux(1:nSubsample_aux:end), mrRate_clu(1:nSubsample_aux:end,iClu), 'k.');
        xlabel(vcLabel_aux);
        ylabel('Firing Rate (Hz)');
        grid on;
    end %for
end %func
