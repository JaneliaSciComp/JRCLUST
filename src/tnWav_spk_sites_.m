%--------------------------------------------------------------------------
% 171201 JJJ: Unique sites handling for diagonal plotting
function tnWav_spk1 = tnWav_spk_sites_(viSpk1, viSites1, S0, fWav_raw_show)
    % reorder tnWav1 to viSites1
    % P = get0_('P');
    % if nargin<3, fWav_raw_show = P.fWav_raw_show; end
    if nargin<3, S0 = []; end
    if isempty(S0), S0 = get(0, 'UserData'); end
    if nargin<4, fWav_raw_show = get_set_(S0.P, 'fWav_raw_show', 0); end

    % unique exception handling %171201 JJJ
    [viSites1_uniq, ~, viiSites1_uniq] = unique(viSites1);
    if numel(viSites1_uniq) ~= numel(viSites1)
        tnWav_spk11 = tnWav_spk_sites_(viSpk1, viSites1_uniq, S0, fWav_raw_show);
        tnWav_spk1 = tnWav_spk11(:,viiSites1_uniq,:);
        return;
    end

    [viSite_spk, P] = deal(S0.viSite_spk, S0.P);
    tnWav = get_spkwav_(P, fWav_raw_show);
    nT_spk = size(tnWav, 1);
    nSpk1 = numel(viSpk1);
    viSites_spk1 = viSite_spk(viSpk1);
    viSites_spk_unique = unique(viSites_spk1);
    tnWav_spk1 = zeros([nT_spk, numel(viSites1), nSpk1], 'like', tnWav);
    for iSite1 = 1:numel(viSites_spk_unique) %only care about the first site
        iSite11 = viSites_spk_unique(iSite1); %center sites group
        viSpk11 = find(viSites_spk1 == iSite11); %dangerous error
        viSites11 = P.miSites(:, iSite11);
        [vlA11, viiB11] = ismember(viSites11, viSites1);
        tnWav_spk1(:,viiB11(vlA11),viSpk11) = tnWav(:,vlA11,viSpk1(viSpk11));
    end
end %func
