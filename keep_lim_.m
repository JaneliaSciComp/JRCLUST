%--------------------------------------------------------------------------
function vi = keep_lim_(vi, lim)
    vi = vi(vi>=lim(1) & vi <= lim(end));
end
