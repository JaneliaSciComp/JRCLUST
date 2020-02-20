function obj = calculateTemplateSim(obj,flagged,max_lag_s)
    if nargin<3
        max_lag_s=5e-4; % 0.5 ms
        if nargin<2
            flagged=[];
        end
    end
    max_lag_samples = round(max_lag_s * obj.hCfg.sampleRate);
    if ~isempty(flagged)
        obj.templateSim(flagged,flagged) = calculate_simscore(obj.meanWfLocalRaw(:,:,flagged),max_lag_samples,obj.hCfg.siteNeighbors(:,obj.clusterSites(flagged)));
    else
        obj.templateSim = calculate_simscore(obj.meanWfLocalRaw,max_lag_samples,obj.hCfg.siteNeighbors(:,obj.clusterSites));
    end
end