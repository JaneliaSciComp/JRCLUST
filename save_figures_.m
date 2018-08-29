%--------------------------------------------------------------------------
function save_figures_(vcExt)
    % bottom to top left to right
    global fDebug_ui

    vcPrefix = sprintf('jrc3_%s_', datestr(now, 'yymmdd-HHMM'));
    csAns = inputdlg_('Figure name prefix', 'Save figure set', 1, {vcPrefix});
    if isempty(csAns), return; end
    vcPrefix = csAns{1};

    figTags = get0_('figTags');
    fprintf('Saving figures...\n'); t1=tic;
    for iFig=1:numel(figTags)
        %     hFig1 = P.vhFig(iField);
        hFig1 = getCachedFig(figTags{iFig});
        if ~ishandle(hFig1), continue; end
        vcFile1 = [vcPrefix, get(hFig1, 'Tag'), vcExt];
        if fDebug_ui==1, continue; end
        %     if getOr([], 'fDebug_ui', 0), continue; end %skip saving for debugging
        switch lower(vcExt)
            case '.fig'
            savefig(hFig1, vcFile1, 'compact');
            otherwise
            saveas(hFig1, vcFile1, vcExt(2:end));
        end
        fprintf('\t%s\n', vcFile1);
    end
    fprintf('\ttook %0.1fs\n', toc(t1));
end % function
