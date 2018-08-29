%--------------------------------------------------------------------------
function P = upgrade_param_(S0, P0)
    % Upgrade parameter struct P based on the old value P0

    % upgrade raw
    P = S0.P;
    nSamples_raw = S0.traceDims(1);
    nSamples_raw_P = diff(S0.P.spkLim_raw) + 1;
    if nSamples_raw_P ~= nSamples_raw
        spkLim_raw = P.spkLim * 2;
        nSamples_raw_P = diff(spkLim_raw) + 1;
        if nSamples_raw_P == nSamples_raw
            P.spkLim_raw = spkLim_raw;
            P.spkLim_raw_factor = 2;
            P.spkLim_raw_ms = [];
        elseif nSamples_raw == S0.waveformDims(1)
            P.spkLim_raw = P.spkLim;
            P.spkLim_raw_factor = 1;
            P.spkLim_raw_ms = [];
            P.spkLim_factor_merge = 1;
        else
            P.spkLim_raw = P0.spkLim_raw;
            P.spkLim_raw_factor = P.spkLim_raw(1) /  P.spkLim(1);
            P.spkLim_raw_ms = [];
            %         disperr_('Dimension of spkLim_raw is inconsistent between S0 and P');
        end
    end

    % upgrade P.maxSite and P.nSites_ref (v3.1.1 to v3.1.9)
    nSites_spk = S0.traceDims(2);
    nSites_spk_P = P.maxSite * 2 + 1;
    if nSites_spk_P ~= nSites_spk
        P.maxSite = (nSites_spk-1)/2;
        P.nSites_ref = nSites_spk - S0.featureDims(1);
        P.miSites = findNearSites_(P.mrSiteXY, P.maxSite, P.viSiteZero, P.viShank_site);
    end
end % function
