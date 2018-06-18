%--------------------------------------------------------------------------
function [viTime1, viSite1, viSpk1] = trimSpikes_(nlim_show)
    error('not implemented');

    viTime1=[];  viSite1=[];
    if ~isfield(P, 'viSpk') || ~isfield(P, 'viSite'), return; end
    ilim = round(P.tlim * P.sRateHz);
    viSpk1 = find(P.viSpk >= ilim(1) & P.viSpk < ilim(end));
    viTime1 = P.viSpk(viSpk1) - ilim(1);
    viSite1 = P.viSite(viSpk1);
end %func
