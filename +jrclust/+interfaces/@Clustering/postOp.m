function postOp(obj, updateMe)
    %POSTOP Call this after a merge or split operation
    % update counts, center sites, remove empty clusters
    obj.refresh(1, updateMe);

    if ~isempty(updateMe)
        % update cluster waveforms and distance
        obj.updateWaveforms(updateMe);

        % update cluster positions
        obj.computeCentroids(updateMe);
    
        % compute quality scores for new clusters
        obj.computeQualityScores(updateMe);
    end
end