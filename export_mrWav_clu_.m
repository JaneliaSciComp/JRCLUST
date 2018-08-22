%--------------------------------------------------------------------------
function export_mrWav_clu_(h,e)
    % Export selected cluster waveforms (iCopy, iPaste) set in the GUI

    S0 = get(0, 'UserData');
    P = S0.P;
    S_clu = S0.S_clu;
    viSite1 = P.miSites(:,S_clu.clusterSites(S0.primarySelectedCluster));
    mrWav_clu1 = S_clu.tmrWav_clu(:,viSite1,S0.primarySelectedCluster);
    eval(sprintf('mrWav_clu%d = mrWav_clu1;', S0.primarySelectedCluster));
    eval(sprintf('assignWorkspace_(mrWav_clu%d);', S0.primarySelectedCluster));
    if ~isempty(S0.secondarySelectedCluster);
        mrWav_clu2 = S_clu.tmrWav_clu(:,viSite1,S0.secondarySelectedCluster);
        eval(sprintf('mrWav_clu%d = mrWav_clu2;', S0.secondarySelectedCluster));
        eval(sprintf('assignWorkspace_(mrWav_clu%d);', S0.secondarySelectedCluster));
    end
end %func
