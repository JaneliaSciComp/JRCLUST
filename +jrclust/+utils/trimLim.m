function lim1 = trimLim(lim1, lim0)
    %TRIMLIM Find intersection of two limit ranges
    dx = diff(lim1);

    if lim1(1) < lim0(1)
        lim1 = lim0(1) + [0, dx];
    end
    if lim1(2) > lim0(2)
        lim1 = lim0(2) + [-dx, 0];
    end
    lim1(1) = max(lim1(1), lim0(1));
    lim1(2) = min(lim1(2), lim0(2));
end
