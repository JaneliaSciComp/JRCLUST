%--------------------------------------------------------------------------
% 9/27/17 JJJ: Created and tested
function vcFilter = get_filter_(P)
    if isfield(P, 'vcFilter')
        vcFilter = P.vcFilter;
    else
        if get_(P, 'nDiff_filt') > 0
            vcFilter = 'sgdiff'; %sgdiff?
        else
            vcFilter = 'bandpass';
        end
    end
end %func
