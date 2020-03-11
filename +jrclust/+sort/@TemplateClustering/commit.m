function success = commit(obj, spikeClusters, metadata, msg)
    %COMMIT Commit a modification of clustering to history log
    success = commit@jrclust.interfaces.Clustering(obj, spikeClusters, metadata, msg);
    if success
        if isempty(fieldnames(metadata)) % initial commit                     
            obj = obj.calculateSimScore([]);          
        elseif isfield(metadata, 'unitCount') && any(isnan(metadata.unitCount))
            flagged = isnan(metadata.unitCount);        
            obj = obj.calculateSimScore(flagged); 
        end
    end
end

