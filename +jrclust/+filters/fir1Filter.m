function samplesIn = fir1Filter(samplesIn, order, passband)
    %FIR1FILTER Convolve samples with a FIR bandpass filter of given order
    ff = single(fir1(order, passband));
    for i = 1:size(samplesIn, 2)
        samplesIn(:, i) = conv(samplesIn(:, i), ff, 'same');
    end
end

