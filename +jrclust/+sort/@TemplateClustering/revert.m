function success = revert(obj, revertTo)
    %REVERT Summary of this function goes here
    success = revert@jrclust.interfaces.Clustering(obj, revertTo);

    if success
        obj.templatesByCluster = arrayfun(@(iC) unique(obj.spikeTemplates(obj.spikesByCluster{iC})), ...
                                  1:obj.nClusters, 'UniformOutput', 0);

        templateSim = zeros(obj.nClusters);
        for iCluster = 1:obj.nClusters
            iTemplates = obj.templatesByCluster{iCluster}; % unique template indices for spikes in this cluster

            % compute cluster sim score, Phy style
            sims = max(obj.sRes.simScore(iTemplates, :), [], 1);

            for jCluster = iCluster:obj.nClusters
                jTemplates = obj.templatesByCluster{jCluster};
                templateSim(iCluster, jCluster) = max(sims(jTemplates));
                templateSim(jCluster, iCluster) = templateSim(iCluster, jCluster);
            end
        end

        obj.templateSim = templateSim;
    end
end
