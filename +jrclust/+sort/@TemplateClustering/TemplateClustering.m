classdef TemplateClustering < jrclust.interfaces.Clustering
    %TEMPLATECLUSTERING A Kilosort(2) clustering    
    %% TEMPLATE-MATCHING PROPERTIES
    properties (Dependent, SetObservable)
        amplitudes;         % extracted amplitude for each spike
        simScore;           % template-template correlation matrix
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
            dpFields = jsondecode(fread(fid, inf, '*char')');
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
        function vals = get.amplitudes(obj)
            if isfield(obj.sRes, 'amplitudes')
                vals = obj.sRes.amplitudes;
            else
                vals = [];
            end
        end
        function set.amplitudes(obj, vals)
            obj.sRes.amplitudes = vals;
        end

        function vals = get.simScore(obj)
            if isfield(obj.sRes, 'simScore')
                vals = obj.sRes.simScore;
            else
                vals = [];
            end
        end
        function set.simScore(obj, vals)
            obj.sRes.simScore = vals;
        end

        function vals = get.spikeTemplates(obj)
            if isfield(obj.sRes, 'spikeTemplates')
                vals = obj.sRes.spikeTemplates;
            else
                vals = [];
            end
        end
        function set.spikeTemplates(obj, vals)
            obj.sRes.spikeTemplates = vals;
        end
    end
end

