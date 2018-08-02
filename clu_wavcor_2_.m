%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function vrWavCor2 = clu_wavcor_2_(ctmrWav_clu, clusterSites, iClu2, cell_5args)
    [P, vlClu_update, mrWavCor0, cviShift1, cviShift2] = deal(cell_5args{:});
    nClu = numel(clusterSites);
    iSite_clu2 = clusterSites(iClu2);
    if iSite_clu2==0 || isnan(iSite_clu2), vrWavCor2 = []; return; end
    viSite2 = P.miSites(:,iSite_clu2);
    maxDist_site_um = get_set_(P, 'maxDist_site_merge_um', 35);
    viClu1 = find(ismember(clusterSites, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));

    vrWavCor2 = zeros(nClu, 1, 'single');
    viClu1(viClu1 >= iClu2) = []; % symmetric matrix comparison
    if isempty(viClu1), return; end
    cmrWav_clu2 = cellfun(@(x)x(:,viSite2,iClu2), ctmrWav_clu, 'UniformOutput', 0);
    ctmrWav_clu1 = cellfun(@(x)x(:,viSite2,viClu1), ctmrWav_clu, 'UniformOutput', 0);
    for iClu11 = 1:numel(viClu1)
        iClu1 = viClu1(iClu11);
        if ~vlClu_update(iClu1) && ~vlClu_update(iClu2)
            vrWavCor2(iClu1) = mrWavCor0(iClu1, iClu2);
        else
            iSite_clu1 = clusterSites(iClu1);
            if iSite_clu1==0 || isnan(iSite_clu1), continue; end
            if iSite_clu1 == iSite_clu2
                cmrWav_clu2_ = cmrWav_clu2;
                cmrWav_clu1_ = cellfun(@(x)x(:,:,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            else
                viSite1 = P.miSites(:,iSite_clu1);
                viSite12 = find(ismember(viSite2, viSite1));
                if isempty(viSite12), continue; end
                cmrWav_clu2_ = cellfun(@(x)x(:,viSite12), cmrWav_clu2, 'UniformOutput', 0);
                cmrWav_clu1_ = cellfun(@(x)x(:,viSite12,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
            end
            vrWavCor2(iClu1) = maxCor_drift_2_(cmrWav_clu2_, cmrWav_clu1_, cviShift1, cviShift2);
        end
    end %iClu2 loop
end %func
