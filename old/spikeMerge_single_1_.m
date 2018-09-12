%--------------------------------------------------------------------------
function [viSpkA, vrSpkA, viSiteA] = spikeMerge_single_1_(viSpk, vrSpk, viSite, iSite, P)

    maxDist_site_um = get_set_(P, 'maxDist_site_um', 50);
    nPad_pre = get_set_(P, 'nPad_pre', 0);
    nRefrac = int32(abs(P.spkRefrac));
    spkLim = [-nRefrac, nRefrac];

    % filter local sites only
    vii1 = find(viSite == iSite);
    viSpk1 = viSpk(vii1);
    vrSpk1 = vrSpk(vii1);
    viSiteNear = findNearSite_(P.mrSiteXY, iSite, maxDist_site_um);
    vi2 = find(ismember(viSite, viSiteNear));
    viSpk2 = viSpk(vi2);
    vrSpk2 = vrSpk(vi2);
    viSite2 = viSite(vi2);
    n2 = numel(viSpk2);

    vlSpk1 = false(size(viSpk1));
    i2prev = 1;
    for iSpk1=1:numel(viSpk1)
        iSpk11 = viSpk1(iSpk1);
        rSpk11 = vrSpk1(iSpk1);

        % check for duplicate detection. search nearby spikes
        spkLim11 = iSpk11 + spkLim;
        [vii2, i2prev] = findRange_(viSpk2, spkLim11(1), spkLim11(2), i2prev, n2);
        if numel(vii2)==1, vlSpk1(iSpk1) = 1; continue; end %no other spikes detected
        vrSpk22 = vrSpk2(vii2);

        %invalid if larger (more negative) spike found
        if any(vrSpk22 < rSpk11), continue; end %wouldn't work for biopolar spike

        % check for equal amplitude, pick first occured
        vii22Eq = find(vrSpk22 == rSpk11);
        if numel(vii22Eq) > 1
            viSpk22 = viSpk2(vii2);
            viSpk222 = viSpk22(vii22Eq);
            minTime = min(viSpk222);
            if minTime ~= iSpk11, continue; end
            if sum(minTime == viSpk222) > 1 %pick lower site
                viSite22 = viSite2(vii2);
                if min(viSite22(vii22Eq)) ~= iSite, continue; end
            end
        end
        vlSpk1(iSpk1) = 1; % set this spike as valid
    end %for

    % Trim
    viiSpk1 = find(vlSpk1); %speed up since used multiple times
    viSpkA = viSpk1(viiSpk1);
    vrSpkA = vrSpk1(viiSpk1);
    viSiteA = repmat(int32(iSite), size(viiSpk1));
    % fprintf('.'); %progress. how about within parfor?

    % apply spike merging on the same site
    % if ~isempty(refrac_factor) && refrac_factor~=0
    %     nRefrac2 = int32(double(nRefrac) * refrac_factor);
    %     [viSpkA, vrSpkA, viSiteA] = spike_refrac_(viSpkA, vrSpkA, viSiteA, nRefrac2); %same site spikes
    % end
end %func
