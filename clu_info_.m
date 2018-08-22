%--------------------------------------------------------------------------
function clu_info_(S0)
    % This also plots cluster position
    if nargin<1, S0 = get(0, 'UserData'); end
    P = S0.P; S_clu = S0.S_clu;
    mh_info = get_tag_('mh_info', 'uimenu');
    S_clu1 = get_cluInfo_(S0.primarySelectedCluster);
    if ~isempty(S0.secondarySelectedCluster)
        S_clu2 = get_cluInfo_(S0.secondarySelectedCluster);
        vcLabel = sprintf('Unit %d "%s" vs. Unit %d "%s"', ...
        S0.primarySelectedCluster, S_clu.clusterNotes{S0.primarySelectedCluster}, ...
        S0.secondarySelectedCluster, S_clu.clusterNotes{S0.secondarySelectedCluster});
        set(mh_info, 'Label', vcLabel);
    else
        S_clu2 = [];
        vcLabel = sprintf('Unit %d "%s"', S0.primarySelectedCluster, S_clu.clusterNotes{S0.primarySelectedCluster});
        set(mh_info, 'Label', vcLabel);
    end
    plot_FigPos_(S_clu1, S_clu2);
end %func
