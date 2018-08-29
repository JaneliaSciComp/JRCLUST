%--------------------------------------------------------------------------
function vc = if_on_off_(vc, cs)
    if ischar(cs), cs = {cs}; end
    vc = ifeq_(ismember(vc, cs), 'on', 'off');
end % function
