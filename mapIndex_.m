%--------------------------------------------------------------------------
function [vi, nClu, viA] = mapIndex_(vi, viA, viB)
    % change the index of vi according to the map (viA)

    if nargin<2, viA = setdiff(unique(vi), 0); end %excl zero
    if nargin<3, viB = 1:numel(viA); end
    nClu = viB(end);
    viAB(viA) = viB; %create a translation table A->B
    vl = vi>0;
    vi(vl) = viAB(vi(vl)); %do not map zeros
end % function
