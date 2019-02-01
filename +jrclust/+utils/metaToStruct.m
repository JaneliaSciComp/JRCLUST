function S = metaToStruct(filename)
    %METATOSTRUCT read a SpikeGLX meta file and convert to struct

    fid = fopen(filename, 'r');
    keysvals = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', 0);
    fclose(fid);

    keys = keysvals{1};
    vals = keysvals{2};
    S = struct();

    for i = 1:numel(keys)
        key = keys{i};
        if key(1) == '~'
            key(1) = [];
        end
        try
            eval(sprintf('%s = ''%s'';', key, vals{i}));
            eval(sprintf('num = str2double(%s);', key));
            if ~isnan(num)
                eval(sprintf('%s = num;', key));
            end
            eval(sprintf('S = setfield(S, ''%s'', %s);', key, key));
        catch ME
            error('error reading %s: %s', filename, ME.message);
        end
    end
end

