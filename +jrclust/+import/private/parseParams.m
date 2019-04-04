function params = parseParams(paramFile)
    %PARSEPARAMS Parse a Phy params.py file, returning a struct
    params = [];
    if exist(paramFile, 'file') ~= 2
        return;
    end

    params = struct();
    keysVals = cellfun(@(line) strsplit(line, '='), jrclust.utils.readLines(paramFile), 'UniformOutput', 0);
    for i = 1:numel(keysVals)
        kv = cellfun(@strip, keysVals{i}, 'UniformOutput', 0);
        key = kv{1}; val = kv{2};

        switch key
            case 'dat_path' % get full path to dat_path
                val = jrclust.utils.absPath(strip(val, ''''), fileparts(paramFile));

            case {'n_channels_dat', 'offset', 'sample_rate'}
                val = str2double(val);

            case 'hp_filtered'
                val = strcmp(val, 'True');
                
            case 'dtype'
                val = strip(val, '''');
        end

        params.(key) = val;
    end
end

