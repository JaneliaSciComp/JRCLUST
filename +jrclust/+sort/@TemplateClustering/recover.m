function flag = recover(obj, ask)
    %RECOVER Inconsistent data recovery, with a hammer
    if nargin < 2
        ask = 0;
    end

    flag = recover@jrclust.interfaces.Clustering(obj, ask);
    if flag < 1
        return;
    end

    inconsistentFields = obj.inconsistentFields();

    if flag == 2 || ismember('templatesByCluster', inconsistentFields)
        obj.templatesByCluster = arrayfun(@(iC) unique(obj.spikeTemplates(obj.spikesByCluster{iC})), ...
                                  1:obj.nClusters, 'UniformOutput', 0);
    end

    if flag == 2 || ismember('templateSim', inconsistentFields)
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

