%--------------------------------------------------------------------------
function exit_manual_(src, event)
    try
        if ~ishandle(src), return; end
        if ~isvalid(src), return; end
        S0 = get(0, 'UserData');
        P = S0.P;
        %     if ~get_set_([], 'fDebug_ui', 0)
        fExit = save_manual_(P);
        if ~fExit, return; end
        if ~isfield(S0, 'csFig')
            S0.csFig = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigWavCor', 'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
        end
        deleteMany(get_fig_all_(S0.csFig), src);
        tryClose(get_fig_('FigTrial'));
        tryClose(get_fig_('FigTrial_b'));
        tryClose(get_fig_('FigAux'));
    catch
        disperr_();
        close(src);
    end
    set(0, 'UserData', []); % clear previous
end %func
