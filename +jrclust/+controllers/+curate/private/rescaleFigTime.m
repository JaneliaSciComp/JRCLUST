function rescaleFigTime(hFigTime, timeScale)
    %set_fig_maxAmp_('FigTime', timeScale);
    hFigTime.axSet('YLim', [0, 1]*timeScale);
    imrect_set_(hFigTime, 'hRect', [], [0, timeScale]);

    % switch lower(P.vcFet_show)
    %     case {'vpp', 'vmin'} %voltage feature
    %         vcYlabel = sprintf('Site %d (\\mu%s)', iSite, P.vcFet_show);
    %     otherwise %other feature options
    %         vcYlabel = sprintf('Site %d (%s)', iSite, P.vcFet_show);
    % end
    % ylabel(S_fig.hAx, vcYlabel);

end %func
