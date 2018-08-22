%--------------------------------------------------------------------------
function auto_scale_proj_time_(S0, fPlot)
    % auto-scale and refresh
    if nargin < 1
        S0 = get(0, 'UserData');
    end

    if nargin < 2
        fPlot = 0;
    end

    autoscale_pct = getOr(S0.P, 'autoscale_pct', 99.5);
    [hFig_proj, S_fig_proj] = getCachedFig('FigProj');
    [mrMin0, mrMax0, mrMin1, mrMax1, mrMin2, mrMax2] = fet2proj_(S0, S_fig_proj.sitesOfInterest);

    if isempty(mrMin2) || isempty(mrMax2)
        cmrAmp = {mrMin1, mrMax1};
    else
        cmrAmp = {mrMin1, mrMax1, mrMin2, mrMax2};
    end
    S_fig_proj.maxAmp = max(cellfun(@(x) quantile(x(:), autoscale_pct/100), cmrAmp));
    set(hFig_proj, 'UserData', S_fig_proj);

    % Update time
    [hFig_time, S_fig_time] = getCachedFig('FigTime');
    iSite = S0.S_clu.clusterSites(S0.primarySelectedCluster);

    [vrFet1, vrTime1, vcYlabel, viSpk1] = getFet_site_(iSite, S0.primarySelectedCluster, S0); % plot primarySelectedCluster
    if isempty(S0.iCluPaste)
        cvrFet = {vrFet1};
    else
        [vrFet2, vrTime2, vcYlabel, viSpk2] = getFet_site_(iSite, S0.iCluPaste, S0); % plot primarySelectedCluster
        cvrFet = {vrFet1, vrFet2};
    end

    S_fig_time.maxAmp = max(cellfun(@(x) quantile(x(:), autoscale_pct/100), cvrFet));
    set(hFig_time, 'UserData', S_fig_time);

    % plot
    if fPlot
        keyPressFcn_cell_(getCachedFig('FigWav'), {'j', 't'}, S0);
    else
        rescale_FigProj_(S_fig_proj.maxAmp, hFig_proj, S_fig_proj, S0);
        rescale_FigTime_(S_fig_time.maxAmp, S0, S0.P);
    end
end %func
