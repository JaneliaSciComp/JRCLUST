%--------------------------------------------------------------------------
function [S_clu, iClu2] = clu_reorder_(S_clu, iClu1, iClu2)
    % Move iClu2 next to iClu1
    % from end to location. iClu1: place next to iClu1
    % reorder clu location
    % cluster is already appended
    if nargin<2,  iClu2 = S_clu.nClu; end

    if nargin < 2
        iClu1 = find(S_clu.viSite_clu < S_clu.viSite_clu(end), 1, 'last');
        if isempty(iClu1), return; end
    end
    iClu2 = iClu1+1; %move one right to iClu1
    if iClu2 == S_clu.nClu, return; end %if no change in position return

    vlAdd = S_clu.viClu>iClu1 & S_clu.viClu < S_clu.nClu;
    S_clu.viClu(S_clu.viClu == S_clu.nClu) = iClu2;
    S_clu.viClu(vlAdd) = S_clu.viClu(vlAdd) + 1;

    viMap_clu = 1:S_clu.nClu;
    viMap_clu(iClu2) = S_clu.nClu;
    viMap_clu(iClu1+2:end) = (iClu1+1):(S_clu.nClu-1);

    S_clu = S_clu_select_(S_clu, viMap_clu);
end %func
