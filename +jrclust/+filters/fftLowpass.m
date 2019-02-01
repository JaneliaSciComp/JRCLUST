function samplesOut = fftLowpass(samplesIn, fCut, sampleRate)
    %FFTLOWPASS Low-pass filter samplesIn
    if isempty(fCut)
        samplesOut = samplesIn;
        return;
    end

    shape = size(samplesIn);

    samplesX = fft(reshape(samplesIn, shape(1), []));
    freqs_ = fftshift((1:size(samplesX, 1))/size(samplesX, 1)) - .5;

    keepMe = abs(freqs_) <= fCut/sampleRate/2;
    samplesX(~keepMe, :) = 0;

    samplesOut = real(reshape(ifft(samplesX), shape));
end
