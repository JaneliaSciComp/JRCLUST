%--------------------------------------------------------------------------
function vc = if_on_off_(vc, cs)
    if ischar(cs), cs = {cs}; end
    vc = jrclust.utils.ifEq(ismember(vc, cs), 'on', 'off');
end %func
