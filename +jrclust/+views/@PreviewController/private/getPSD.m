function [psdPower, psdFreq] = getPSD(traces, sampleRate, nSkip)
    %GETPSD Compute power spectral density of traces
    if nSkip > 1
        traces = traces(1:nSkip:end, :);
    end

    n = size(traces, 1);
    n1 = round(n/2);

    psdPower = fft(jrclust.utils.meanSubtract(single(traces)));
    psdPower = jrclust.utils.pow2db(abs(psdPower(2:n1 + 1, :))) / n;

    if nargout > 1
        psdFreq = sampleRate*(1:n1)'/n;
    end
end