%--------------------------------------------------------------------------
function tr1 = fft_lowpass_(tr, fc, sampleRateHz)
    if isempty(fc), tr1 = tr; return; end
    dimm = size(tr);

    mr = fft(reshape(tr, dimm(1), []));
    vrFreq = fftshift((1:size(mr,1))/size(mr,1)) - .5;
    vlUse = abs(vrFreq) <= fc/sampleRateHz/2;
    mr(~vlUse,:) = 0;
    tr1 = real(reshape(ifft(mr), dimm));
end %func
