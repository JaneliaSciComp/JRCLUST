%--------------------------------------------------------------------------
function S_clu = split_clu_(iClu1, vlIn)
    % split cluster.
    figure_wait_(1); drawnow;
    [P, S_clu, viSite_spk] = get0_('P', 'S_clu', 'viSite_spk');
    hMsg = msgbox_open_('Splitting...');
    figure(getCachedFig('FigWav'));

    % create a new cluster (add at the end)
    n2 = sum(vlIn); %number of clusters to split off
    iClu2 = max(S_clu.viClu) + 1;

    % update cluster count and index
    S_clu.nClusters = double(iClu2);
    S_clu.vnSpk_clu(iClu1) = S_clu.vnSpk_clu(iClu1) - n2;
    S_clu.vnSpk_clu(iClu2) = sum(vlIn);
    viSpk1 = find(S_clu.viClu==iClu1);
    viSpk2 = viSpk1(vlIn);
    viSpk1 = viSpk1(~vlIn);
    [iSite1, iSite2] = deal(mode(viSite_spk(viSpk1)), mode(viSite_spk(viSpk2)));
    if iSite1 > iSite2 % order by the cluster site location
        [iSite2, iSite1] = deal(iSite1, iSite2);
        [viSpk2, viSpk1] = deal(viSpk1, viSpk2);
        vlIn = ~vlIn;
    end
    [S_clu.cviSpk_clu{iClu1}, S_clu.cviSpk_clu{iClu2}] = deal(viSpk1, viSpk2);
    [S_clu.viSite_clu(iClu1), S_clu.viSite_clu(iClu2)] = deal(iSite1, iSite2);

    try % erase annotaiton
        S_clu.csNote_clu{iClu1} = '';
        S_clu.csNote_clu{end+1} = ''; %add another entry
    catch
    end
    S_clu.viClu(viSpk2) = iClu2; %change cluster number
    S_clu = S_clu_update_(S_clu, [iClu1, iClu2], P);

    % Bring the new cluster right next to the old one using index swap
    [S_clu, iClu2] = clu_reorder_(S_clu, iClu1);

    % update all the other views
    [S_clu, S0] = S_clu_commit_(S_clu, 'split_clu_');
    plot_FigWav_(S0); %redraw plot
    plotFigWavCor(S0);
    save_log_(sprintf('split %d', iClu1), S0); %@TODO: specify which cut to use

    % select two clusters being split
    button_CluWav_simulate_(iClu1, iClu2);
    tryClose(hMsg);
    fprintf('%s [W] splited Clu %d\n', datestr(now, 'HH:MM:SS'), iClu1);
    figure_wait_(0);
end %func
