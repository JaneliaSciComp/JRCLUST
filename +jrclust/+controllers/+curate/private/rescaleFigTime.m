function rescaleFigTime(hFigTime, timeScale)
    %set_fig_maxAmp_('FigTime', timeScale);
    YLim = hFigTime.axApply('default', @get, 'YLim');
    hFigTime.axApply('default', @set, 'YLim', [0, YLim(2)*timeScale]);
    imrect_set_(hFigTime, 'hRect', [], [0, YLim(2)*timeScale]);

    % switch lower(P.vcFet_show)
    %     case {'vpp', 'vmin'} %voltage feature
    %         vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.vcFet_show);
    %     otherwise %other feature options
    %         vcYlabel = sprintf('Site %d (%s)', iSite, P.vcFet_show);
    % end
    % ylabel(S_fig.hAx, vcYlabel);

end
