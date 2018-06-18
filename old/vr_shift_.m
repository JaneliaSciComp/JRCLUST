%--------------------------------------------------------------------------
function mr = vr_shift_(vr, viShift);
    % viShift = [0, -1,-.5,.5,1]; %[0, -.5, .5]
    % viShift = [0, -1:.25:-.25,.25:.25:1];
    nShifts = numel(viShift);
    vr = vr(:);
    vi0 = (1:numel(vr))';
    mr = zeros(numel(vr), nShifts, 'like', vr);
    for iShift = 1:nShifts
        dn = viShift(iShift);
        if dn~=0
            mr(:,iShift) = interp1(vi0, vr, vi0+dn, 'pchip', 'extrap');
        else
            mr(:,iShift) = vr;
        end
    end
end %func
