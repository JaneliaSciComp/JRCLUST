%--------------------------------------------------------------------------
function traces = tnWav2uV_(traces, P)
    % Scale traces and subtract off their means.

    if nargin < 2
        P = get0_('P');
    end

    traces = bit2uV_(traces, P);
    traces = meanSubtract(traces);
end % function
