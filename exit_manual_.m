%--------------------------------------------------------------------------
function exit_manual_(src, event)
    try
        if ~ishandle(src), return; end
        if ~isvalid(src), return; end
        S0 = get(0, 'UserData');
        P = S0.P;
        fExit = save_manual_(P);
        if ~fExit, return; end
        if ~isfield(S0, 'figTags')
            S0.figTags = {'FigPos', 'FigMap', 'FigTime', 'FigWav', 'FigClusterCor', ...
                'FigProj', 'FigRD', 'FigCorr', 'FigIsi', 'FigHist'};
        end
        deleteMany(figuresByTag(S0.figTags), src);
        tryClose(figuresByTag('FigTrial'));
        tryClose(figuresByTag('FigTrial_b'));
        tryClose(figuresByTag('FigAux'));
    catch
        disperr_();
        close(src);
    end
    set(0, 'UserData', []); % clear previous
end %func
