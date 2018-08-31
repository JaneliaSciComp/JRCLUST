%--------------------------------------------------------------------------
function  [viTop, vrTop] = find_topn_(vr, nMax, vi)
    if nargin<3, vi = 1:numel(vr); end
    nMax = min(nMax, numel(vi));
    [~, viSrt] = sort(vr(vi), 'descend');
    viTop = vi(viSrt(1:nMax));
    vrTop= vr(viTop);
end % function
