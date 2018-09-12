%--------------------------------------------------------------------------
function restore_clu_(varargin)
    % restore last deleted. most negative clu is last deleted
    % error('to be fixed. Fix centroid code');
    S0 = get(0, 'UserData');
    [P, S_clu] = deal(S0.P, S0.S_clu);
    iClu_del = min(S_clu.viClu);
    if iClu_del >=0, msgbox_('Deleted cluster is not found'); return; end

    figure_wait_(1);
    % if deleted add a clu at the end and zoom at it
    % change clusters
    iClu_new = double(max(S_clu.viClu) + 1);
    S_clu.viClu(S_clu.viClu == iClu_del) = iClu_new;
    S_clu.nClu = iClu_new;
    S_clu = S_clu_update_(S_clu, iClu_new, P);
    S_clu.csNote_clu{end+1} = '';
    [S_clu, iClu_new] = clu_reorder_(S_clu);

    % update all the other views
    % delete_multi_(S0.vhPlot, S0.vhText);
    % S0.S_clu = S_clu; set(0, 'UserData', S0);
    [S_clu, S0] = S_clu_commit_(S_clu);
    plot_FigWav_(S0); %redraw plot
    plot_FigWavCor_(S0);
    % S0 = set0_(mrWavCor);
    set(0, 'UserData', S0);

    % append to the end for now
    button_CluWav_simulate_(iClu_new);
    keyPressFcn_cell_(get_fig_cache_('FigWav'), 'z');
    fprintf('%s [W] Restored Clu %d\n', datestr(now, 'HH:MM:SS'), iClu_new);
    figure_wait_(0);
end
