function params = getDefaultParams(merge)
    %GETDEFAULTPARAMS Get default JRCLUST params
    %   merge - bool,optional: merge common and advanced parameters if true
    if nargin < 1
        merge = 1;
    end

    fid = fopen(fullfile(jrclust.utils.basedir(), 'json', 'params.json'), 'r', 'n', 'UTF-8');
    fstr = fread(fid, '*char')';
    fclose(fid);

    params = jsondecode(fstr);
    if merge
        params = jrclust.utils.mergeStructs(params.commonParameters, params.advancedParameters);
    end
end

