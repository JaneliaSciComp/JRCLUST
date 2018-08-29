%--------------------------------------------------------------------------
function vc = str_toggle_(vc, vc1, vc2)
    % toggle vc1 to vc2
    if strcmpi(vc, vc1)
        vc = vc2;
    else
        vc = vc1;
    end
end % function
