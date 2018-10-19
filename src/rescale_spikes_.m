%--------------------------------------------------------------------------
function rescale_spikes_(hSpkAll, maxAmp_prev, P)
    S = get(hSpkAll, 'UserData');
    [~, S_fig] = get_fig_cache_('FigWav');
    S0 = get(0, 'UserData');
    cvrY = S.cvrY;
    cviSite = S.cviSite;
    % nSamples = diff(P.spkLim)+1;
    scale = S_fig.maxAmp / maxAmp_prev;
    for iClu=1:numel(cvrY)
        viSite1 = cviSite{iClu};
        nSites1 = numel(viSite1);
        trY = reshape(cvrY{iClu}, [], nSites1, S.vnSpk(iClu));
        for iSite1 = 1:nSites1
            y_off = viSite1(iSite1);
            trY(:,iSite1,:) = (trY(:,iSite1,:) - y_off) / scale + y_off;
        end
        cvrY{iClu} = trY(:);
    end
    S.cvrY = cvrY;
    set(hSpkAll, 'YData', cell2mat_(cvrY), 'UserData', S);
end
