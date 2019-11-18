function updateSubset(obj, showSubset, silent)
%UPDATESUBSET Update the subset of units to show
if nargin < 3
    silent = 0;
end

showSubset = intersect(1:obj.hClust.nClusters, showSubset);

if isempty(showSubset) && ~silent
    jrclust.utils.qMsgBox('Subset is empty; selecting all units');
    showSubset = obj.showSubset;
elseif jrclust.utils.isEqual(1:obj.hClust.nClusters, showSubset(:)') && ~silent
    jrclust.utils.qMsgBox('All units selected');
end

if ~jrclust.utils.isEqual(obj.showSubset(:), showSubset(:))
    obj.showSubset = showSubset;

    obj.replot();
end
end % function

