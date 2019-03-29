function refresh(obj, doRemoveEmpty, updateMe)
    %REFRESH Summary of this function goes here
    if nargin < 2
        doRemoveEmpty = 0;
    end
    if nargin < 3
        updateMe = [];
    end

    refresh@jrclust.interfaces.Clustering(obj, doRemoveEmpty, updateMe);

    obj.templatesByCluster = arrayfun(@(iCluster) unique(obj.spikeTemplates(obj.spikesByCluster{iCluster})), ...
                                      1:obj.nClusters, 'UniformOutput', 0);

    if ~isempty(obj.nClusters)
        templateSim_ = zeros(obj.nClusters);
        for iCluster = 1:obj.nClusters
            iTemplates = obj.templatesByCluster{iCluster}; % unique template indices for spikes in this cluster

            % compute cluster sim score, Phy style
            sims = max(obj.sRes.simScore(iTemplates, :), [], 1);

            for jCluster = iCluster:obj.nClusters
                jTemplates = obj.templatesByCluster{jCluster};
                templateSim_(iCluster, jCluster) = max(sims(jTemplates));
                templateSim_(jCluster, iCluster) = templateSim_(iCluster, jCluster);
            end
        end

        obj.templateSim = templateSim_;
    end
end

