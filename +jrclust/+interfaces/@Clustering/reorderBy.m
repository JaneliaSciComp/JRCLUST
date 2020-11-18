function [success, argsort] = reorderBy(obj, by)
%REORDERBY Arrange cluster ID numbers by some criterion.
%   Currently supported criteria are:
%   - Y + X: sort by centroid in (Y, X) lexicographical order
%   - clusterSites: sort by site assignment
%   - any other property which is sortable
if nargin < 2 || isempty(by) || (~strcmp(by, 'Y + X') && ~isprop(obj, by))
    by = 'clusterSites';
end

argsort = 1:obj.nClusters;
if strcmpi(by, 'Y + X') && ~isempty(obj.clusterCentroids)
    [~, argsort] = sort(sum(obj.clusterCentroids, 2), 'ascend');
elseif isprop(obj, by)
    try
        [~, argsort] = sort(obj.(by), 'ascend');
    catch ME
        warning('Failed to reorder: %s', ME.message);
        success = 0;
        return;
    end
end

if issorted(argsort)
    success = 1;
    return;
end

success = obj.reorderMultiple(1:obj.nClusters, argsort);
end %fun