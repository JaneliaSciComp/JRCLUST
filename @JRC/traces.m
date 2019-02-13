function traces(obj)
    %TRACES Show traces
    if isempty(obj.hCfg)
        error('Please specify a config file');
    end

    if isempty(obj.res)
        obj.loadFiles();
    end

    hTraces = jrclust.views.TracesController(obj.hCfg);
    if numel(obj.args) > 1
        recID = str2double(obj.args{2});
        if isnan(recID)
            recID = [];
        end
    else
        recID = [];
    end

    hTraces.show(recID, 0, obj.hClust);
end