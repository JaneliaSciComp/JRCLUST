%--------------------------------------------------------------------------
function vr = linmap_(vr, lim1, lim2, fSat)
    if nargin< 4
        fSat = 0;
    end
    if numel(lim1) == 1, lim1 = [-abs(lim1), abs(lim1)]; end
    if numel(lim2) == 1, lim2 = [-abs(lim2), abs(lim2)]; end

    if fSat
        vr(vr>lim1(2)) = lim1(2);
        vr(vr<lim1(1)) = lim1(1);
    end
    if lim1(1)==lim1(2)
        vr = vr / lim1(1);
    else
        vr = interp1(lim1, lim2, vr, 'linear', 'extrap');
    end
end %func
