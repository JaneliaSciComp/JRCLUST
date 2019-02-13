function saveFigures(obj, ext)
    %SAVEFIGURES Save figures
    prefix = sprintf('jrc4_%s_', datestr(now, 'yyyy-mm-dd-HHMM'));
    dlgAns = inputdlg('Figure name prefix', 'Save figure set', 1, {prefix});
    if isempty(dlgAns)
        return;
    end
    prefix = dlgAns{1};

    if obj.hCfg.verbose
        fprintf('Saving figures...\n');
        t1 = tic;
    end

    figKeys = keys(obj.hFigs);
    for iFig = 1:numel(figKeys)
        hFig = obj.hFigs(figKeys{iFig});
        if ~hFig.isReady
            continue;
        end

        filename = fullfile(obj.hCfg.outputDir, [prefix, hFig.figApply(@get, 'Tag'), ext]);
        if strcmpi(ext, '.fig')
            savefig(hFig.hFig, filename, 'compact');
        else
            saveas(hFig.hFig, filename, ext(2:end));
        end

        if obj.hCfg.verbose
            fprintf('\t%s\n', filename);
        end
    end

    if obj.hCfg.verbose
        fprintf('\ttook %0.1fs\n', toc(t1));
    end
end