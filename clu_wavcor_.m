%--------------------------------------------------------------------------
% 10/27/17 JJJ: distance-based neighboring unit selection
function vrWavCor2 = clu_wavcor_(ctmrWav_clu, cviSite_clu, P, cell_5args, iClu2)

    [vlClu_update, cviShift1, cviShift2, mrWavCor0, fMode_cor] = deal(cell_5args{:});
    if numel(cviSite_clu) == 1
        viSite_clu = cviSite_clu{1};
        fUsePeak2 = 0;
    else
        [viSite_clu, viSite2_clu, viSite3_clu] = deal(cviSite_clu{:});
        fUsePeak2 = 1;
    end
    nClu = numel(viSite_clu);
    iSite_clu2 = viSite_clu(iClu2);
    if iSite_clu2==0 || isnan(iSite_clu2), vrWavCor2 = []; return; end
    viSite2 = P.miSites(:,iSite_clu2);
    % if fMaxSite_excl, viSite2 = viSite2(2:end); end
    if fUsePeak2
        viClu1 = find(viSite_clu == iSite_clu2 | viSite2_clu == iSite_clu2 | viSite3_clu == iSite_clu2 | ...
        viSite_clu == viSite2_clu(iClu2) | viSite_clu == viSite3_clu(iClu2)); %viSite2_clu == viSite2_clu(iClu2)
    else
        %     maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
        maxDist_site_um = get_set_(P, 'maxDist_site_merge_um', 35);
        viClu1 = find(ismember(viSite_clu, findNearSite_(P.mrSiteXY, iSite_clu2, maxDist_site_um)));
    end

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
            iSite_clu1 = viSite_clu(iClu1);
            if iSite_clu1==0 || isnan(iSite_clu1), continue; end
            if iSite_clu1 == iSite_clu2
                cmrWav_clu2_ = cmrWav_clu2;
                cmrWav_clu1_ = cellfun(@(x)x(:,:,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
                diff=false; % TW; appears to have no effect -- acl
            else
                viSite1 = P.miSites(:,iSite_clu1);
                viSite12 = find(ismember(viSite2, viSite1));
                viSite23 = find(ismember(viSite1, viSite2)); % TW; appears to have no effect -- acl
                if isempty(viSite12), continue; end
                cmrWav_clu2_ = cellfun(@(x)x(:,viSite12), cmrWav_clu2, 'UniformOutput', 0);
                cmrWav_clu1_ = cellfun(@(x)x(:,viSite12,iClu11), ctmrWav_clu1, 'UniformOutput', 0);
                diff=true; % TW; appears to have no effect -- acl
            end
            vrWavCor2(iClu1) = maxCor_drift_(cmrWav_clu2_, cmrWav_clu1_, cviShift1, cviShift2, fMode_cor);
        end
    end %iClu2 loop
end %func
