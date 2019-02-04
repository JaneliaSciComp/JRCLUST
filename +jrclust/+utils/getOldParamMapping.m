function [old2new, new2old] = getOldParamMapping()
    %GETOLDPARAMMAPPING Get the old-style parameter mapping
    fid = fopen(fullfile(jrclust.utils.basedir, 'params.json'), 'r', 'n', 'UTF-8');
    fstr = fread(fid, '*char')';
    fclose(fid);

    params = jsondecode(fstr);
    old2new = params.old2new;

    if nargout > 1
        new2old = struct();
        oldFn = fieldnames(old2new);
        for i = 1:numel(oldFn)
            new2old.(old2new.(oldFn{i})) = oldFn{i};
        end
    end
end

