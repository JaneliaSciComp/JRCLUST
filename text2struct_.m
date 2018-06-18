%--------------------------------------------------------------------------
% 8/2/17 JJJ: Documentation and test
function S = text2struct_(vcFname)
    % Convert text file to struct

    fid = fopen(vcFname, 'r');
    mcFileMeta = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
    fclose(fid);
    csName = mcFileMeta{1};
    csValue = mcFileMeta{2};
    S = struct();
    for i=1:numel(csName)
        vcName1 = csName{i};
        if vcName1(1) == '~', vcName1(1) = []; end
        try
            eval(sprintf('%s = ''%s'';', vcName1, csValue{i}));
            eval(sprintf('num = str2double(%s);', vcName1));
            if ~isnan(num)
                eval(sprintf('%s = num;', vcName1));
            end
            eval(sprintf('S = setfield(S, ''%s'', %s);', vcName1, vcName1));
        catch
            fprintf('%s = %s error\n', csName{i}, csValue{i});
        end
    end
end %func
