function resetPositions(obj)
    obj.figApply(@(hFig) hFig.resetPos());
    obj.allToForeground();
end