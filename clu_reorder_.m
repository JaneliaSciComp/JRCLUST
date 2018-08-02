%--------------------------------------------------------------------------
function [S_clu, iClu2] = clu_reorder_(S_clu, iClu1, iClu2)
    % Move iClu2 next to iClu1
    % from end to location. iClu1: place next to iClu1
    % reorder clu location
    % cluster is already appended
    if nargin<2,  iClu2 = S_clu.nClusters; end

    if nargin < 2
        iClu1 = find(S_clu.clusterSites < S_clu.clusterSites(end), 1, 'last');
        if isempty(iClu1), return; end
    end
    iClu2 = iClu1+1; %move one right to iClu1
    if iClu2 == S_clu.nClusters, return; end %if no change in position return

    vlAdd = S_clu.spikeClusters>iClu1 & S_clu.spikeClusters < S_clu.nClusters;
    S_clu.spikeClusters(S_clu.spikeClusters == S_clu.nClusters) = iClu2;
    S_clu.spikeClusters(vlAdd) = S_clu.spikeClusters(vlAdd) + 1;

    viMap_clu = 1:S_clu.nClusters;
    viMap_clu(iClu2) = S_clu.nClusters;
    viMap_clu(iClu1+2:end) = (iClu1+1):(S_clu.nClusters-1);

    S_clu = S_clu_select_(S_clu, viMap_clu);
end %func
