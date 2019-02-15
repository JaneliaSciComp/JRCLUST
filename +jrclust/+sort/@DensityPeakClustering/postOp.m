function postOp(obj, updateMe)
    %POSTOP Extend default Clustering postOp method
    if ~isempty(updateMe)
        arrayfun(@obj.rmRefracSpikes, updateMe);
    end

    postOp@jrclust.interfaces.Clustering(obj, updateMe);
end