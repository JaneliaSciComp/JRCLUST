%--------------------------------------------------------------------------
function [mnWav2, cviSite_mean] = meanSite_drift_(mnWav1, P, viSite_repair)
    % [mnWav1, cviSite_mean] = meanSite_drift_(mnWav1, P)
    % [mnWav1, cviSite_mean] = meanSite_drift_(mnWav1, P, viSite_repair) %bad site repair

    % mnWav1 = sites_repair_(mnWav1, P); % must repair site to perform vertical averaging
    % this corrects for viSiteZero automatically
    nSites = size(mnWav1,2);
    nCols = nColumns_probe_(P);
    viSiteZero = get_(P, 'viSiteZero');
    cviSite_mean = cell(1, nSites);
    viSites = 1:nSites;
    fSingleShank = isSingleShank_(P);
    if nargin<3
        viSite_repair = 1:nSites;
    else
        viSite_repair = toRow_(viSite_repair);
    end
    mnWav2 = zeros(size(mnWav1), 'like', mnWav1);
    for iSite = viSite_repair
        if fSingleShank
            viSite_shank1 = viSites; %faster
        else
            viSite_shank1 = viSites(P.viShank_site == P.viShank_site(iSite));
        end
        for iNeigh=1:4
            viSite_mean1 = iSite + [-1,0,0,1] * nCols * iNeigh;
            viSite_mean1 = viSite_mean1(ismember(viSite_mean1, viSite_shank1));
            viSite_mean1(ismember(viSite_mean1, viSiteZero)) = [];
            if ~isempty(viSite_mean1), break; end
        end
        cviSite_mean{iSite} = viSite_mean1;
        mnWav2(:, iSite) = mean(mnWav1(:, viSite_mean1), 2);
    end
end %func

%% local functions
function flag = isSingleShank_(P)
    viShank_site = get_(P, 'viShank_site');
    if isempty(viShank_site)
        flag = 1;
    else
        flag = numel(unique(viShank_site)) == 1;
    end
end

function nCols = nColumns_probe_(P)
    % Checkerboard four-column is considered as two column probe since
    % two sites per vertical step
    viShank_site = get_(P, 'viShank_site');
    if ~isempty(viShank_site)
        viSites = find(P.viShank_site == P.viShank_site(1));
        vrSiteY = P.mrSiteXY(viSites,2);
    else
        vrSiteY = P.mrSiteXY(:,2);
    end
    vrSiteY_unique = unique(vrSiteY);
    vnSites_group = hist(vrSiteY, vrSiteY_unique);
    nCols = median(vnSites_group);
end %func
