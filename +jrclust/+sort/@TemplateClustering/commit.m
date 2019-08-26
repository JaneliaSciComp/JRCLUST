function success = commit(obj, spikeClusters, metadata, msg)
    %COMMIT Commit a modification of clustering to history log
    success = commit@jrclust.interfaces.Clustering(obj, spikeClusters, metadata, msg);

    if success
        if 1==1 % initial commit -- AGB forcing it to enter this if statement fixes the problem of template similarity plot not updating
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
        elseif isfield(metadata, 'unitCount') && any(isnan(metadata.unitCount))
            flagged = find(isnan(metadata.unitCount));
            obj.templatesByCluster(flagged) = arrayfun(@(iC) unique(obj.spikeTemplates(obj.spikesByCluster{iC})), ...
                flagged, 'UniformOutput', 0);

            % update cluster templates
            for i = 1:numel(flagged)
                iCluster = flagged(i);
                iTemplates = obj.templatesByCluster{iCluster}; % unique template indices for spikes in this cluster

                % compute cluster sim score, Phy style
                sims = max(obj.sRes.simScore(iTemplates, :), [], 1);

                for jCluster = 1:obj.nClusters
                    jTemplates = obj.templatesByCluster{jCluster};
                    obj.templateSim(iCluster, jCluster) = max(sims(jTemplates));
                    obj.templateSim(jCluster, iCluster) = obj.templateSim(iCluster, jCluster);
                end
            end
        end
    end
end

