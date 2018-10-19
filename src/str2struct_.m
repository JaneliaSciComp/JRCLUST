%--------------------------------------------------------------------------
% 8/14/17 JJJ: Created
function [S, vcErr] = str2struct_(vc)
    % Numerical expressions only for now

    vc = vc';
    vc = vc(:)';
    vcErr = '';
    cc_ = textscan(strtrim(vc), '%s%s', 'Delimiter', '=;',  'ReturnOnError', false);
    csName = cc_{1};
    csValue = cc_{2};

    S = struct();
    for i=1:numel(csName)
        vcName_ = strtrim(csName{i});
        if vcName_(1) == '~', vcName_(1) = []; end
        try
            eval(sprintf('%s = %s;', vcName_, csValue{i}));
            eval(sprintf('num = str2double(%s);', vcName_));
            if ~isnan(num)
                eval(sprintf('%s = num;', vcName_));
            end
            eval(sprintf('S = setfield(S, ''%s'', %s);', vcName_, vcName_));
        catch
            vcErr = lasterr();
            S = [];
            return;
            %         fprintf('%s = %s error\n', csName{i}, csValue{i});
        end
    end

end %func
