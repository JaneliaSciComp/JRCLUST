function obj = calculateSimScore(obj,flagged,max_lag_s)
    if nargin<3
        max_lag_s=5e-4; % 0.5 ms
        if nargin<2
            flagged=[];
        end
    elseif strcmp(obj.hCfg.simMode,'template')
        warning('max_lag_s not used for template-based similarity scores.');
    end
    switch obj.hCfg.simMode
        case 'waveform'
            % recalculate similarity scores between clusters using max
            % cross-correlation (up to max_lag_s) between mean waveforms
            max_lag_samples = round(max_lag_s * obj.hCfg.sampleRate);
            if ~isempty(flagged)
                obj.templateSim(flagged,flagged) = jrclust.sort.TemplateClustering.waveformSimScore(...
                    obj.meanWfLocalRaw(:,:,flagged),max_lag_samples,obj.hCfg.siteNeighbors(:,obj.clusterSites(flagged)));
            else
                obj.templateSim = jrclust.sort.TemplateClustering.waveformSimScore(...
                    obj.meanWfLocalRaw,max_lag_samples,obj.hCfg.siteNeighbors(:,obj.clusterSites));
            end
        case 'template'
            % use precomputed template-based similarity scores         
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