%--------------------------------------------------------------------------
function clu_info_(S0)
    % This also plots cluster position
    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    hMenuInfo = get_tag_('hMenuInfo', 'uimenu');
    S_clu1 = get_cluInfo_(S0.iCluCopy);
    if ~isempty(S0.iCluPaste)
        S_clu2 = get_cluInfo_(S0.iCluPaste);
        vcLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', ...
        S0.iCluCopy, S_clu.csNote_clu{S0.iCluCopy}, ...
        S0.iCluPaste, S_clu.csNote_clu{S0.iCluPaste});
        set(hMenuInfo, 'Label', vcLabel);
    else
        S_clu2 = [];
        vcLabel = sprintf('Unit %d "%s"', S0.iCluCopy, S_clu.csNote_clu{S0.iCluCopy});
        set(hMenuInfo, 'Label', vcLabel);
    end
    plot_FigPos_(S_clu1, S_clu2);
end %func
