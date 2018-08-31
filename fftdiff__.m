%--------------------------------------------------------------------------
function mnWav1 = fftdiff__(mnWav, freqLim_)
    % apply fft to diffrentiate
    % mnWav = gather_(mnWav);

    n = size(mnWav,1);
    % n1 = round(n/2*freqLim_(1));
    % n2 = round(n/2*diff(freqLim_));

    n1 = round(n/2 * freqLim_(2));
    npow2 = 2^nextpow2(n);
    % w = single([linspace(0, 1, n2), linspace(1, 0, n2)])';
    % w = [zeros(n1, 1, 'single'); w; zeros(npow2-2*n1-4*n2, 1, 'single'); -w; zeros(n1, 1, 'single')];
    % w = single(pi*1i) * w;
    w = single(pi*1i) * single([linspace(0, 1, n1), linspace(1, -1, npow2-2*n1), linspace(-1, 0, n1)]');
    mnWav1 = real(ifft(bsxfun(@times, fft(single(mnWav), npow2), w), 'symmetric'));
    mnWav1 = cast(mnWav1(1:n,:), class_(mnWav));
end % function
