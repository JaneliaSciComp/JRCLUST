function allToForeground(obj)
    %ALLTOFOREGROUND Bring all figures to the foreground
    obj.figApply(@(hFig) hFig.toForeground());
end