%--------------------------------------------------------------------------
function vals = linmap(vals, oldLimits, newLimits, saturate)
    % maintain relative distribution of vals while mapping into newLimits
    if nargin < 4
        saturate = 0;
    end

    if numel(oldLimits) == 1
        oldLimits = [-abs(oldLimits), abs(oldLimits)];
    end
    if numel(newLimits) == 1
        newLimits = [-abs(newLimits), abs(newLimits)];
    end

    if saturate
        vals(vals > oldLimits(2)) = oldLimits(2);
        vals(vals < oldLimits(1)) = oldLimits(1);
    end

    if oldLimits(1) == oldLimits(2) % ignore newLimits and just rescale
        vals = vals/oldLimits(1);
    else
        vals = interp1(oldLimits, newLimits, vals, 'linear', 'extrap');
    end
end %func
