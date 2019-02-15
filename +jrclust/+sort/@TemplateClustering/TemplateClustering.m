classdef TemplateClustering < jrclust.interfaces.Clustering
    %TEMPLATECLUSTERING A Kilosort(2) clustering
    
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

            obj.clearNotes();
            obj.refresh(1, []);
            commitMsg = sprintf('%s;initial import', datestr(now, 31));
            obj.commit(commitMsg);
        end
    end

    %% GETTERS/SETTERS
    
end

