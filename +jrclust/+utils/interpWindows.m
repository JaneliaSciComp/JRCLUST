function [spikeWindows, indices] = interpWindows(spikeWindows, nIterp)
    %INTERPWINDOWS Interpolate windows subpixel-wise using FFT method
    if nIterp <= 1
        return;
    end

    nSamples = size(spikeWindows, 1);
    % interpolation done on first non-singleton dimension but we assume
    % columnwise; force with 3rd arg
    spikeWindows = interpft(spikeWindows, size(spikeWindows, 1)*nIterp, 1);

    indices = 1:(1/nIterp):nSamples;
    if ndims(spikeWindows) == 3
        spikeWindows = spikeWindows(1:numel(indices), :, :);
    else
        spikeWindows = spikeWindows(1:numel(indices), :);
    end
end
