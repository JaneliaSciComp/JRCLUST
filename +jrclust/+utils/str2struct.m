function S = str2struct(strval)
    %STR2STRUCT Convert a string to a struct
    %   (only numerical expressions supported for now)

    strval = strval';
    strval = strval(:)';

    keysVals = textscan(strtrim(strval), '%s%s', 'Delimiter', '=;',  'ReturnOnError', 0);
    structKeys = keysVals{1};
    structVals = keysVals{2};

    S = struct();
    for i = 1:numel(structKeys)
        fn = strtrim(structKeys{i});

        if fn(1) == '~'
            fn(1) = [];
        end

        try
            eval(sprintf('%s = %s;', fn, structVals{i}));
            eval(sprintf('num = str2double(%s);', fn));

            if ~isnan(num)
                eval(sprintf('%s = num;', fn));
            end

            eval(sprintf('S = setfield(S, ''%s'', %s);', fn, fn));
        catch ME
            warning('str2struct failed: %s', ME.message);
            S = [];
            return;
        end
    end

end %func
