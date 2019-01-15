function tc = trueClass(vals)
    %TRUECLASS Return class of vals, or, in the case of gpuArrays, the
    %   underlying class
    if isempty(vals)
        tc = class(jrclust.utils.tryGather(vals));
    else
        tc = class(jrclust.utils.tryGather(vals(1)));
    end
end
