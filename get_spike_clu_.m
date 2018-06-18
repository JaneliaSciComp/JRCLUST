%--------------------------------------------------------------------------
function [vrX_spk, vrY_spk, vrA_spk, viClu_spk, vrT_spk] = get_spike_clu_(S_clu, viSites, viClu_plot)

    if nargin<3, viClu_plot=[]; end
    S0 = get(0, 'UserData');
    P = S0.P;

    nSpikes = numel(S0.viTime_spk);
    [vrX_spk, vrY_spk, vrA_spk, viClu_spk, vrT_spk] = deal(nan(nSpikes, 1, 'single'));
    for iSite = viSites
        viSite1 = P.miSites(:, iSite);
        viSpk1 = find(S0.viSite_spk == iSite);
        viClu1 = S_clu.viClu(viSpk1);
        if ~isempty(viClu_plot)
            vl1_clu = ismember(viClu1, viClu_plot);
            viSpk1 = viSpk1(vl1_clu);
            viClu1 = viClu1(vl1_clu);
        end
        if isempty(viSpk1), continue; end
        viTime1 = S0.viTime_spk(viSpk1);
        %     [vrX_spk(viSpk1), vrY_spk(viSpk1), vrA_spk(viSpk1)] = ...
        %         spikePos_(viSpk1, viSite1, P);
        viClu_spk(viSpk1) = viClu1;
        vrT_spk(viSpk1) = viTime1;
    end %for

    % select
    vl_spk = ~isnan(vrY_spk);
    vrA_spk = (vrA_spk(vl_spk));
    vrY_spk = vrY_spk(vl_spk) / P.um_per_pix;
    vrX_spk = vrX_spk(vl_spk) / P.um_per_pix;
    vrT_spk = single(vrT_spk(vl_spk)) / P.sRateHz;
    viClu_spk = viClu_spk(vl_spk);

    % sort by amplitude
    [vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk] = ...
    sort_ascend_(vrA_spk, vrX_spk, vrY_spk, vrT_spk, viClu_spk);
end %func
