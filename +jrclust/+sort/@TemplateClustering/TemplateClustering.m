classdef TemplateClustering < jrclust.interfaces.Clustering
    %TEMPLATECLUSTERING A Kilosort(2) clustering    
    %% TEMPLATE-MATCHING PROPERTIES
    properties (Dependent, SetObservable)
        spikeTemplates;     % template assignments for each spike
    end

    properties (SetObservable)
        templatesByCluster; % cell array of unique template values by cluster
        templateSim;        % template-based similarity score
    end

    %% LIFECYCLE
    methods
        function obj = TemplateClustering(sRes, dRes, hCfg)
            %TEMPLATECLUSTERING
            fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'TemplateClustering.json'), 'r');
            dpFields = jsondecode(fread(fid, inf, '*char'));
            fclose(fid);
            obj.unitFields.vectorFields = [obj.unitFields.vectorFields; dpFields.vectorFields];

            obj.sRes = sRes;
            obj.dRes = dRes;
            obj.hCfg = hCfg;

            obj.spikeClusters = sRes.spikeClusters;

            obj.clearNotes();
            obj.refresh(1, []);
            commitMsg = sprintf('%s;initial import', datestr(now, 31));
            obj.commit(commitMsg);
        end
    end

    %% UTILITY METHODS
    methods (Access=protected, Hidden)
        nMerged = mergeBySim(obj);
    end

    %% GETTERS/SETTERS
    methods
        function vals = get.spikeTemplates(obj)
            if ~isempty(obj.sRes)
                vals = obj.sRes.spikeTemplates;
            else
                vals = [];
            end
        end
    end
end

