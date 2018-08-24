%--------------------------------------------------------------------------
function traces = bit2uV_(traces, P)
    % Scale traces by normalized uV_per_bit factor

    if nargin < 2
        P = get0_('P');
    end

    switch lower(get_filter_(P))
        case 'sgdiff'
            norm = sum((1:P.nDiff_filt).^2) * 2;

        case 'ndiff'
            norm = 2;

        otherwise
            norm = 1;
    end

    traces = single(traces)*single(P.uV_per_bit/norm);
end
