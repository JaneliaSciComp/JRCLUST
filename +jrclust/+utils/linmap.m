function vals = linmap(vals, oldLim, newLim)
    %LINMAP Rescale vals occurring within oldLim to newLim, saturating at
    %the boundaries
    if numel(oldLim) == 1
        oldLim = abs(oldLim)*[-1, 1];
    end

    % saturate at the boundaries
    vals(vals > oldLim(2)) = oldLim(2);
    vals(vals < oldLim(1)) = oldLim(1);

    if all(oldLim == 0) % nothing to rescale, all is 0 now
        return;
    end

    if oldLim(1) == oldLim(2) % ignore newLim and just rescale 
        vals = vals / oldLim(1);
    else
        vals = interp1(oldLim, newLim, vals, 'linear', 'extrap');
    end
end