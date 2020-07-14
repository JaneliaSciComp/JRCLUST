function success = revert2(obj, toIndex)
%REVERT2 Revert some number of entries in the history.
%   For every delete, do an undelete. For every merge, do a split, for
%   every split, a merge. Each undo operation gets committed to the
%   history.
% toIndex doesn't refer to an entry we can revert to
if toIndex < 0 || toIndex > obj.nEdits
    error('revert2 called with a bad index %d', toIndex);
end


end

