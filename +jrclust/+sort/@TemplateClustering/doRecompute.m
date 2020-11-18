function success = doRecompute(obj)
%DORECOMPUTE Recompute derivative properties of the spike table after an
%edit.
%   After basic doRecompute, update templatesByCluster and template
%   similarity score, then clear obj.recompute.
success = doRecompute@jrclust.interfaces.Clustering(obj);

recompute = obj.recompute(:);
if numel(obj.recompute) == 0
    recompute = (1:obj.nClusters)';
end

% recompute templatesByCluster
if success
    try
        templatesByCluster = cellfun(@(sC) unique(obj.spikeTemplates(sC)), ...
            obj.spikesByCluster(recompute), 'UniformOutput', 0);
        obj.templatesByCluster(recompute) = templatesByCluster;
    catch ME
        warning(ME.message);
        success = 0;
    end
end

% recompute simScore
if success
    templateSim = obj.templateSim;

    try
        for i = 1:numel(recompute)
            iCluster = recompute(i);
            % unique template indices for spikes in this cluster
            iTemplates = obj.templatesByCluster{iCluster};

            % compute cluster sim score, Phy style
            sims = max(obj.sRes.simScore(iTemplates, :), [], 1);

            for jCluster = 1:obj.nClusters
                jTemplates = obj.templatesByCluster{jCluster};
                templateSim(iCluster, jCluster) = max(sims(jTemplates));
                templateSim(jCluster, iCluster) = templateSim(iCluster, jCluster);
            end
        end

        obj.templateSim = templateSim;
    catch ME
        warning(ME.message);
        success = 0;
    end
end

if success
    obj.recompute = [];
end
end %fun

