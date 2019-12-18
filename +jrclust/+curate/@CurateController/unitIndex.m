function unitIndices = unitIndex(obj, unitIDs)
%UNITINDEX Get index into showSubset of unitIDs specified
unitIndices = find(ismember(obj.showSubset, unitIDs));
end

