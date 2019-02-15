function revert(obj, revertTo)
    %REVERT Delete history
    if revertTo < 0 || revertTo >= obj.nEdits
        return;
    end

    obj.editSeek(revertTo);
    obj.history(revertTo+2:end, :) = [];
    obj.removeEmptyClusters();
    obj.computeCentroids();
    obj.computeQualityScores([]);
    obj.clearNotes();
    obj.editPos = revertTo;
end